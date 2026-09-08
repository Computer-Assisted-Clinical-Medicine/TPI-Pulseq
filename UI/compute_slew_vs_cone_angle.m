function diag = compute_slew_vs_cone_angle(theta_used_rad, G_amp, k0, kmax, p, ...
    gamma, sys, ramp_mode, ramp_stretch, transition_mode, transition_stretch)
% ==========================================================================
% COMPUTE_SLEW_VS_CONE_ANGLE  Full 1-deg sweep of theoretical peak slew
%                              rate across all cone angles, with markers
%                              for the angles actually compiled.
%
% Inputs:
%   theta_used_rad - polar angles [rad] actually compiled (e.g.
%                    theta_cones(cone_indices) after equator filtering);
%                    used only to flag matching sweep points.
%   G_amp,k0,kmax,p,gamma,sys,ramp_mode,ramp_stretch,transition_mode,
%   transition_stretch - same design inputs as design_tpi_projection.
%
% Output diag fields:
%   theta_deg              - dense pole-safe sweep grid [deg] (0.5:1:179.5)
%   Sx, Sy, Sz              - per-axis peak slew across the sweep [T/m/s]
%   Snorm                   - vector-norm peak slew [T/m/s]
%   Saxis                   - per-axis max, i.e. the Pulseq criterion [T/m/s]
%   Sanalytical             - Boada Eq. [11] theta-dependent bound
%   limit                   - hardware slew limit [T/m/s]
%   used_mask               - true at sweep points nearest a compiled cone
%   used_theta_deg/used_Snorm/used_Saxis - subsets for marker plotting
% 
% Author:  Valentin Jost
% Date:    2026-09
% ==========================================================================

if nargin < 10 || isempty(transition_mode),    transition_mode = 'none'; end
if nargin < 11 || isempty(transition_stretch), transition_stretch = 1.0; end

theta_deg = (0.5:1:179.5)';
theta_rad = deg2rad(theta_deg);
NS = length(theta_deg);

S_x = zeros(NS,1); S_y = zeros(NS,1); S_z = zeros(NS,1);
S_vec = zeros(NS,1); S_ax = zeros(NS,1);

for n = 1:NS
    [~, ~, ~, wfm_n, ~] = design_tpi_projection( ...
        theta_rad(n), G_amp, k0, kmax, p, gamma, sys, ramp_mode, ramp_stretch, ...
        transition_mode, transition_stretch, true);
    if wfm_n.invalid
        S_x(n) = NaN; S_y(n) = NaN; S_z(n) = NaN; S_vec(n) = NaN; S_ax(n) = NaN;
    else
        S_x(n)   = wfm_n.S_max_x;
        S_y(n)   = wfm_n.S_max_y;
        S_z(n)   = wfm_n.S_max_z;
        S_vec(n) = wfm_n.S_max_constructed;
        S_ax(n)  = wfm_n.S_max_per_axis;
    end
end

S_limit      = sys.maxSlew / gamma;
S_analytical = (gamma * G_amp)^2 ./ (k0 * sin(theta_rad)) / gamma;   % Boada Eq. [11]

used_mask = false(NS,1);
theta_used_deg = rad2deg(theta_used_rad(:));
for u = 1:length(theta_used_deg)
    [~, idx] = min(abs(theta_deg - theta_used_deg(u)));
    used_mask(idx) = true;
end

diag = struct('theta_deg', theta_deg, 'Sx', S_x, 'Sy', S_y, 'Sz', S_z, ...
    'Snorm', S_vec, 'Saxis', S_ax, 'Sanalytical', S_analytical, 'limit', S_limit, ...
    'used_mask', used_mask, ...
    'used_theta_deg', theta_deg(used_mask), ...
    'used_Snorm', S_vec(used_mask), ...
    'used_Saxis', S_ax(used_mask));

end