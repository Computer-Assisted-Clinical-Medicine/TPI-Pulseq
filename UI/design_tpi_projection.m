function [gx, gy, gz, wfm, t_ramp] = design_tpi_projection( ...
    theta_0, G_amp, k0, kmax, p, gamma, sys, ramp_mode, ramp_stretch, ...
    transition_mode, transition_stretch, suppress_warnings)
% ==========================================================================
% DESIGN_TPI_PROJECTION  Build one TPI readout gradient for a single cone.
%
% Constant-polar-angle (cone) trajectory from Boada et al. (1997):
%   Phase 1 - Radial ramp:      k = 0    -> k0    (trapezoid along theta_0)
%   Phase 2 - TPI twist:        k = k0   -> kmax   (Eq. [4], [9])
%   Phase 3 - Ramp-down:        |G| -> 0           (Pulseq block integrity)
%
% Inputs:
%   theta_0            - cone polar angle [rad]
%   G_amp              - constant gradient amplitude [T/m]
%   k0, kmax           - twist onset and max k-space radius [1/m]
%   p                  - twist onset fraction (= k0/kmax)
%   gamma              - gyromagnetic ratio [Hz/T]
%   sys                - Pulseq system limits struct
%   ramp_mode          - 'linear' (default) or 'smooth' (raised-cosine)
%   ramp_stretch       - extra duration multiplier on the |G| ramp (default 1.0)
%   transition_mode    - 'none' (default) or 'taper': raised-cosine ramp
%                        of dphi/dt at the Phase 1/2 boundary, removing the
%                        hard step in angular velocity that otherwise
%                        spikes the local slew rate
%   transition_stretch - multiplier on taper duration (default 1.0)
%   suppress_warnings  - if true, silence slew warnings and return an
%                        'invalid' wfm instead of throwing on hardware
%                        rejection (used for diagnostic angle sweeps)
%
% Outputs:
%   gx, gy, gz - Pulseq mr.makeArbitraryGrad objects ([] if invalid)
%   wfm        - struct with raw waveform, slew metrics, and status flags
%
% Author:  Valentin Jost
% Date:    2026-09
% ==========================================================================

if nargin < 8  || isempty(ramp_mode),          ramp_mode = 'linear';      end
if nargin < 9  || isempty(ramp_stretch),       ramp_stretch = 1.0;       end
if nargin < 10 || isempty(transition_mode),    transition_mode = 'none'; end
if nargin < 11 || isempty(transition_stretch), transition_stretch = 1.0; end
if nargin < 12 || isempty(suppress_warnings),  suppress_warnings = false; end

dt       = sys.gradRasterTime;
G_amp_Hz = gamma * G_amp;   % [Hz/m]

%% ---- Phase 1: radial ramp k = 0 -> k0 ----
n_ramp_steps_linear = ceil(G_amp_Hz / sys.maxSlew / dt);
switch ramp_mode
    case 'linear'
        n_ramp_steps = ceil(n_ramp_steps_linear * ramp_stretch);
        n_ramp       = n_ramp_steps + 1;
        t_ramp       = n_ramp_steps * dt;
        G_p1_ramp    = linspace(0, G_amp, n_ramp);
    case 'smooth'
        n_ramp_steps = ceil(n_ramp_steps_linear * (pi/2) * ramp_stretch);
        n_ramp       = n_ramp_steps + 1;
        t_ramp       = n_ramp_steps * dt;
        tau          = linspace(0, 1, n_ramp);
        G_p1_ramp    = G_amp * (1 - cos(pi * tau)) / 2;
    otherwise
        error('design_tpi_projection: unknown ramp_mode "%s" (use ''linear'' or ''smooth'')', ramp_mode);
end

k_after_ramp = gamma * trapz(G_p1_ramp) * dt;
k_flat1      = max(0, k0 - k_after_ramp);
t_flat1      = ceil(k_flat1 / (gamma * G_amp) / dt) * dt;
n_flat1      = round(t_flat1 / dt);

G_p1      = [G_p1_ramp, G_amp * ones(1, n_flat1)];
k0_actual = gamma * (sum(G_p1) * dt);

Gx_p1 = G_p1 .* sin(theta_0);
Gy_p1 = zeros(1, length(G_p1));
Gz_p1 = G_p1 .* cos(theta_0);
t_phase2_start = length(G_p1) * dt;

%% ---- Phase 2: TPI twist k = k0 -> kmax ----
T_read   = kmax / (gamma * G_amp) * (1 + 2*p^3) / (3*p^2);   % Boada Eq. [14]
T_phase2 = ceil((T_read - t_phase2_start) / dt) * dt;
n_twist  = round(T_phase2 / dt);

t_vec = ((0:n_twist-1) + 0.5) * dt;
k_t   = (3 * gamma * G_amp * k0_actual^2 .* t_vec + k0_actual^3).^(1/3);   % Boada Eq. [4]
dkdt  = gamma * G_amp * k0_actual^2 ./ k_t.^2;

sin_theta_eff = max(sin(theta_0), 1e-3);   % numerical safety near the poles
dphi_dt = sqrt(max(0, (gamma*G_amp)^2 - dkdt.^2)) ./ (k_t .* sin_theta_eff);

% Optional raised-cosine taper of the Phase 1/2 angular-velocity onset,
% independent of ramp_mode/ramp_stretch (which govern only |G|).
switch transition_mode
    case 'none'
        % dphi_dt unchanged: hard onset at the Phase 2 boundary
    case 'taper'
        n_onset_linear = ceil(G_amp_Hz / sys.maxSlew / dt) * 2;
        n_onset = min(n_twist, max(1, ceil(n_onset_linear * transition_stretch)));
        tau_on  = linspace(0, 1, n_onset);
        taper   = ones(1, n_twist);
        taper(1:n_onset) = (1 - cos(pi * tau_on)) / 2;
        dphi_dt = dphi_dt .* taper;
    otherwise
        error('design_tpi_projection: unknown transition_mode "%s" (use ''none'' or ''taper'')', transition_mode);
end

phi_t = cumsum(dphi_dt) * dt;

Gx_twist = (1/gamma) .* (dkdt .* sin(theta_0) .* cos(phi_t) ...
                        - k_t  .* sin(theta_0) .* dphi_dt .* sin(phi_t));
Gy_twist = (1/gamma) .* (dkdt .* sin(theta_0) .* sin(phi_t) ...
                        + k_t  .* sin(theta_0) .* dphi_dt .* cos(phi_t));
Gz_twist = (1/gamma) .* (dkdt .* cos(theta_0));

%% ---- Phase 3: ramp-down to zero ----
Gx_end = Gx_twist(end);
Gy_end = Gy_twist(end);
Gz_end = Gz_twist(end);
G_end_mag_Hz = gamma * sqrt(Gx_end^2 + Gy_end^2 + Gz_end^2);
n_ramp_down_linear = ceil(G_end_mag_Hz / sys.maxSlew / dt);

switch ramp_mode
    case 'linear'
        n_ramp_down = ceil(n_ramp_down_linear * ramp_stretch) + 1;
        ramp_factor = linspace(1, 0, n_ramp_down);
    case 'smooth'
        n_ramp_down = ceil(n_ramp_down_linear * (pi/2) * ramp_stretch) + 1;
        tau_dn      = linspace(0, 1, n_ramp_down);
        ramp_factor = (1 + cos(pi * tau_dn)) / 2;
end

Gx_dn = Gx_end .* ramp_factor;
Gy_dn = Gy_end .* ramp_factor;
Gz_dn = Gz_end .* ramp_factor;

%% ---- Assemble full waveform ----
Gx_wfm = [Gx_p1, Gx_twist, Gx_dn];
Gy_wfm = [Gy_p1, Gy_twist, Gy_dn];
Gz_wfm = [Gz_p1, Gz_twist, Gz_dn];

Gx_Hz = gamma .* Gx_wfm;
Gy_Hz = gamma .* Gy_wfm;
Gz_Hz = gamma .* Gz_wfm;

%% ---- Slew rate metrics (computed before Pulseq object construction) ----
slew_x   = diff(Gx_Hz) / dt;
slew_y   = diff(Gy_Hz) / dt;
slew_z   = diff(Gz_Hz) / dt;
slew_xy  = sqrt(slew_x.^2 + slew_y.^2);
slew_vec = sqrt(slew_x.^2 + slew_y.^2 + slew_z.^2);

S_max_constructed = max(slew_vec) / gamma;                    % vector-norm slew
S_max_ax = max([max(abs(slew_x)), max(abs(slew_y)), max(abs(slew_z))]) / gamma;   % per-axis (Pulseq criterion)

[~, idx_peak] = max(slew_vec);
n_p1 = length(Gx_p1);
n_p2 = n_p1 + length(Gx_twist);
if     idx_peak <= n_p1, peak_phase = 1;
elseif idx_peak <= n_p2, peak_phase = 2;
else,                    peak_phase = 3;
end

if ~suppress_warnings
    if S_max_constructed > sys.maxSlew / gamma
        warning(['design_tpi_projection: vector slew %.1f T/m/s > limit %.1f T/m/s ' ...
                 '(theta=%.1f deg, peak in Phase %d)'], ...
            S_max_constructed, sys.maxSlew/gamma, rad2deg(theta_0), peak_phase);
    end
    if S_max_ax > sys.maxSlew / gamma
        warning(['design_tpi_projection: per-axis slew %.1f T/m/s > limit %.1f T/m/s ' ...
                 '(theta=%.1f deg, Pulseq will reject this waveform)'], ...
            S_max_ax, sys.maxSlew/gamma, rad2deg(theta_0));
    end
end

wfm.S_max_constructed = S_max_constructed;
wfm.S_max_per_axis     = S_max_ax;
wfm.S_max_x            = max(abs(slew_x)) / gamma;
wfm.S_max_y            = max(abs(slew_y)) / gamma;
wfm.S_max_z            = max(abs(slew_z)) / gamma;
wfm.S_max_xy           = max(slew_xy) / gamma;
wfm.slew_profile       = slew_vec / gamma;
wfm.peak_slew_phase    = peak_phase;
wfm.ramp_mode          = ramp_mode;
wfm.t_ramp_up          = t_ramp;
wfm.transition_mode    = transition_mode;
wfm.invalid            = false;

wfm.Gx = Gx_wfm; wfm.Gy = Gy_wfm; wfm.Gz = Gz_wfm;
wfm.t_phase2_start = t_phase2_start;
wfm.phi_t = phi_t;
wfm.k_t   = k_t;

%% ---- Pulseq gradient object construction ----
% mr.makeArbitraryGrad runs its own internal slew check and throws on
% violation. With suppress_warnings=true (diagnostic sweeps over angles
% that may be infeasible) that error is caught instead of propagated.
try
    gx = mr.makeArbitraryGrad('x', Gx_Hz, sys, 'first', 0, 'last', 0);
    gy = mr.makeArbitraryGrad('y', Gy_Hz, sys, 'first', 0, 'last', 0);
    gz = mr.makeArbitraryGrad('z', Gz_Hz, sys, 'first', 0, 'last', 0);
catch ME
    if suppress_warnings
        gx = []; gy = []; gz = [];
        wfm.invalid   = true;
        wfm.error_msg = ME.message;
    else
        rethrow(ME);
    end
end

end