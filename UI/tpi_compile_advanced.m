function out = tpi_compile_advanced(prm)
% ==========================================================================
% TPI_COMPILE_ADVANCED  Compute-only TPI Pulseq compiler (no plotting).
%
% Assembles a full Twisted Projection Imaging sequence for sodium (or
% proton) MRI from the parameter struct prm, following Boada et al.
% (1997). All plotting is deferred to the caller (e.g. TPI_UI).
%
% prm fields (all required):
%   scanner ('Skyra'|'Trio'|'TerraW60'), nuclei ('Na'|'H')
%   fov [m], Nx, TR [s], p, Tread [s], flip_angle_deg
%   ramp_mode ('linear'|'smooth'), ramp_stretch
%   transition_mode ('none'|'taper'), transition_stretch
%   avoid_equator (logical), symmetric_cones (logical)
%   ro_os, TEx (1|2), TE2_delay [s], undersample_factor (>= 1)
%
% Output struct fields:
%   seq, sys, gamma, kmax, k0, G_amp, T_read, S_max_required
%   theta_cones, N_per_cone, N_total
%   gx_test, gy_test, gz_test, wfm_test, theta_test
%   ktraj, t_ktraj, ktraj_adc, t_adc
%   slew_diag (struct: theta/Sx/Sy/Sz/Snorm/Saxis/Sanalytical vectors)
%   log (cell array of status strings for display)
% 
% Author:  Valentin Jost
% Date:    2026-09
% ==========================================================================

log = {};
tic0 = tic;

%% ---- Scanner & system limits ----
switch prm.nuclei
    case 'H';  gamma = 42.576e6;
    case 'Na'; gamma = 11.262e6;
end

safety_margin = 0.95;
switch prm.scanner
    case 'Skyra'
        sys = mr.opts('B0', 3, 'MaxGrad', 45*safety_margin, 'GradUnit', 'mT/m', ...
            'MaxSlew', 200*safety_margin, 'SlewUnit', 'T/m/s', ...
            'rfRingdownTime', 10e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6, 'gamma', gamma);
        rf_dur = 300e-6; sys.rfRingdownTime = 270e-6;
    case 'Trio'
        sys = mr.opts('B0', 3, 'MaxGrad', 40*safety_margin, 'GradUnit', 'mT/m', ...
            'MaxSlew', 200*safety_margin, 'SlewUnit', 'T/m/s', ...
            'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6, 'gamma', gamma);
        rf_dur = 300e-6; sys.rfRingdownTime = 270e-6;
    case 'TerraW60'
        sys = mr.opts('B0', 7, 'MaxGrad', 130*safety_margin, 'GradUnit', 'mT/m', ...
            'MaxSlew', 250*safety_margin, 'SlewUnit', 'T/m/s', ...
            'rfRingdownTime', 10e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6, 'gamma', gamma);
        rf_dur = 500e-6;
end
seq = mr.Sequence(sys);

%% ---- Derived k-space geometry ----
fov = prm.fov; Nx = prm.Nx; TR = prm.TR;
p_twist = prm.p; Tread = prm.Tread; flip_angle_deg = prm.flip_angle_deg;

deltak = 1 / fov;
kmax   = Nx * deltak / 2;
k0     = p_twist * kmax;

ramp_mode           = prm.ramp_mode;
ramp_stretch        = prm.ramp_stretch;
transition_mode     = prm.transition_mode;
transition_stretch  = prm.transition_stretch;
avoid_equator       = prm.avoid_equator;
symmetric_cones     = prm.symmetric_cones;
undersample_factor  = prm.undersample_factor;

ro_os = prm.ro_os;
TEx   = prm.TEx;
TE2_delay = 0; if TEx == 2, TE2_delay = prm.TE2_delay; end

G_amp = kmax / (gamma*Tread) * (1 + 2*p_twist^3) / (3*p_twist^2);   % Boada Eq. [15]
log{end+1} = sprintf('Gradient amplitude:  G = %.2f mT/m (limit %.0f mT/m)', ...
    G_amp*1e3, sys.maxGrad/gamma*1e3);
if gamma*G_amp > sys.maxGrad
    log{end+1} = 'WARNING: G_amp exceeds hardware limit!';
end

T_read = kmax / (gamma*G_amp) * (1 + 2*p_twist^3) / (3*p_twist^2);   % Boada Eq. [14]
log{end+1} = sprintf('Readout time:        T_read = %.3f ms', T_read*1e3);

S_max_required = gamma * G_amp^2 * fov / p_twist;                    % Boada Eq. [13]
log{end+1} = sprintf('Required slew rate:  S = %.1f T/m/s (limit %.0f T/m/s)', ...
    S_max_required, sys.maxSlew/gamma);
if S_max_required > sys.maxSlew/gamma
    log{end+1} = 'WARNING: Slew rate constraint violated. Increase p, reduce G_amp, or use smooth ramp.';
end

%% ---- RF pulse ----
rf = mr.makeBlockPulse(deg2rad(flip_angle_deg), 'Duration', rf_dur, ...
    'system', sys, 'use', 'excitation');

%% ---- Cone distribution (undersampling applied at generation time) ----
[theta_cones, N_per_cone, N_total] = calc_cone_distribution(kmax, fov, p_twist, undersample_factor);

% Fully sampled baseline (undersample_factor = 1), kept for logging only
[theta_cones_full, ~, N_total_full] = calc_cone_distribution(kmax, fov, p_twist, 1);
NR_full = length(theta_cones_full);

%% ---- Symmetric cones: verification/correction pass on N_per_cone ----
if symmetric_cones
    NR_check = length(theta_cones);
    asymmetry_found = false;
    for i = 1:floor(NR_check/2)
        j = NR_check + 1 - i;
        if N_per_cone(i) ~= N_per_cone(j)
            N_sym = max(N_per_cone(i), N_per_cone(j));
            N_per_cone(i) = N_sym;
            N_per_cone(j) = N_sym;
            asymmetry_found = true;
        end
    end
    N_total = sum(N_per_cone);
    if asymmetry_found
        log{end+1} = 'Symmetric cones: rounding asymmetry corrected in mirrored ring pairs.';
    else
        log{end+1} = 'Symmetric cones: ring set already perfectly mirrored (no correction needed).';
    end
end

NR = length(theta_cones);
log{end+1} = sprintf('Number of rings:   NR = %d / %d (undersample factor = %.2f)', ...
    NR, NR_full, undersample_factor);
log{end+1} = sprintf('Total projections: N  = %d / %d', N_total, N_total_full);
log{end+1} = sprintf('Projected DAT:     %.1f / %.1f min (TR=%.0f ms)', ...
    N_total*TR/60, N_total_full*TR/60, TR*1e3);

%% ---- Avoid equator: exact-ring exclusion on the undersampled ring set ----
tol = 1e-6;
dist_to_equator = abs(theta_cones - pi/2);
[min_dist, eq_idx] = min(dist_to_equator);
has_equator_ring = (min_dist < tol);   % true only when NR is odd

exclude_mask = false(NR,1);
if avoid_equator && has_equator_ring
    exclude_mask(eq_idx) = true;
    log{end+1} = sprintf('Equatorial ring detected at theta=%.4f deg (NR=%d, odd) - excluded.', ...
        rad2deg(theta_cones(eq_idx)), NR);
elseif avoid_equator && ~has_equator_ring
    log{end+1} = sprintf('No exact equator ring present (NR=%d, even) - nothing to exclude.', NR);
end

cone_indices = find(~exclude_mask)';
if isempty(cone_indices)
    error('tpi_compile_advanced: no cones remain after exclusion - check avoid_equator settings.');
end
log{end+1} = sprintf('Compiling %d / %d cones (avoid_equator=%d)', length(cone_indices), NR, avoid_equator);

%% ---- Slew diagnostic: full 1-deg sweep, compiled cones marked ----
slew_diag = compute_slew_vs_cone_angle(theta_cones(cone_indices), G_amp, k0, kmax, p_twist, ...
    gamma, sys, ramp_mode, ramp_stretch, transition_mode, transition_stretch);

%% ---- Single reference spoke (example) ----
theta_test = theta_cones(ceil(end/2));
[gx_test, gy_test, gz_test, wfm_test, ~] = design_tpi_projection( ...
    theta_test, G_amp, k0, kmax, p_twist, gamma, sys, ramp_mode, ramp_stretch, ...
    transition_mode, transition_stretch);

log{end+1} = sprintf('Example spoke theta = %.1f deg | S_max constructed = %.1f T/m/s | S_max per-axis = %.1f T/m/s', ...
    rad2deg(theta_test), wfm_test.S_max_constructed, wfm_test.S_max_per_axis);

%% ---- ADC design ----
adc_dur = mr.calcDuration(gx_test);
dwell_nyquist = 1 / (2 * gamma * G_amp * fov);
num_samples   = round(adc_dur / dwell_nyquist * ro_os);
num_samples   = num_samples + mod(num_samples, 2);

adc = mr.makeAdc(num_samples, 'Duration', adc_dur, 'Delay', sys.adcDeadTime, 'system', sys);
block_dur_ref = mr.calcDuration(gx_test, adc);

%% ---- Main assembly loop: undersampled ring set, spokes via mr.rotate ----
TRdelay = TR - (rf_dur + sys.rfRingdownTime + sys.rfDeadTime) - block_dur_ref - TE2_delay;
if TRdelay < 0
    log{end+1} = sprintf('WARNING: TR too short by %.3f ms - clamping delay to 0. Increase TR.', -TRdelay*1e3);
    TRdelay = 0;
end

S_max_all_cones = 0;
worst_theta_deg  = NaN;
worst_phase      = NaN;
worst_metric     = 'none';
n_proj = 0;

for cone_idx = cone_indices
    theta_n  = theta_cones(cone_idx);
    n_spokes = N_per_cone(cone_idx);
    [gx_ref, gy_ref, gz_ref, wfm_ref, ~] = design_tpi_projection( ...
        theta_n, G_amp, k0, kmax, p_twist, gamma, sys, ramp_mode, ramp_stretch, ...
        transition_mode, transition_stretch);

    len_x = length(gx_ref.waveform); len_y = length(gy_ref.waveform); len_z = length(gz_ref.waveform);
    if len_x ~= len_y || len_y ~= len_z || ...
       any(~isfinite(gx_ref.waveform)) || any(~isfinite(gy_ref.waveform)) || any(~isfinite(gz_ref.waveform))
        log{end+1} = sprintf('Skipped degenerate cone at theta=%.2f deg (waveform length/NaN mismatch: %d/%d/%d)', ...
            rad2deg(theta_n), len_x, len_y, len_z);
        continue;
    end

    if wfm_ref.S_max_constructed > S_max_all_cones
        S_max_all_cones = wfm_ref.S_max_constructed;
        worst_theta_deg = rad2deg(theta_n); worst_phase = wfm_ref.peak_slew_phase; worst_metric = 'vector';
    end
    if wfm_ref.S_max_per_axis > S_max_all_cones
        S_max_all_cones = wfm_ref.S_max_per_axis;
        worst_theta_deg = rad2deg(theta_n); worst_phase = wfm_ref.peak_slew_phase; worst_metric = 'per-axis';
    end
    if wfm_ref.S_max_xy > S_max_all_cones
        S_max_all_cones = wfm_ref.S_max_xy;
        worst_theta_deg = rad2deg(theta_n); worst_phase = wfm_ref.peak_slew_phase; worst_metric = 'xy-rotation-invariant';
    end

    if wfm_ref.S_max_xy > sys.maxSlew/gamma
        log{end+1} = sprintf('Skipped cone theta=%.2f deg: xy-rotation-invariant slew %.1f T/m/s exceeds limit %.1f T/m/s (Phase %d)', ...
            rad2deg(theta_n), wfm_ref.S_max_xy, sys.maxSlew/gamma, wfm_ref.peak_slew_phase);
        continue;   % no spoke on this cone would pass
    end

    for j = 1:n_spokes
        n_proj = n_proj + 1;
        phi_j = 2*pi*(j-1)/n_spokes;

        try
            G_rotated = mr.rotate('z', phi_j, gx_ref, gy_ref, gz_ref, 'system', sys);
            seq.addBlock(rf);
            if TEx ~= 1
                seq.addBlock(mr.makeDelay(TE2_delay));
            end
            seq.addBlock(G_rotated{:}, adc);
            seq.addBlock(mr.makeDelay(TRdelay));
        catch ME
            error('tpi_compile_advanced: Pulseq rejected cone theta=%.2f deg, spoke %d/%d (phi=%.1f deg): %s', ...
                rad2deg(theta_n), j, n_spokes, rad2deg(phi_j), ME.message);
        end
    end
end

phase_names = {'Phase 1 (radial ramp)', 'Phase 2 (TPI twist)', 'Phase 3 (ramp-down)'};
phase_str = 'n/a';
if ~isnan(worst_phase), phase_str = phase_names{worst_phase}; end

log{end+1} = sprintf('Compiled %d / %d projections total.', n_proj, N_total_full);
log{end+1} = sprintf('Worst-case slew: %.1f T/m/s (%s metric) at theta=%.1f deg, occurring in %s', ...
    S_max_all_cones, worst_metric, worst_theta_deg, phase_str);

if S_max_all_cones > sys.maxSlew/gamma
    log{end+1} = sprintf('WARNING: worst-case slew exceeds hardware limit (%.0f T/m/s). See cone/phase above.', ...
        sys.maxSlew/gamma);
    if worst_phase == 2
        log{end+1} = 'Note: violation is in Phase 2 (twist) - ramp_mode/ramp_stretch have NO effect here. Increase p, increase Tread, or raise undersample_factor.';
    else
        log{end+1} = 'Note: violation is in Phase 1/3 (ramp) - try ramp_mode=''smooth'' or increase ramp_stretch.';
    end
end

%% ---- k-space trajectory (for plotting) ----
[ktraj_adc, t_adc, ktraj, t_ktraj] = seq.calculateKspacePP();

out = struct('seq', seq, 'sys', sys, 'gamma', gamma, 'kmax', kmax, 'k0', k0, ...
    'G_amp', G_amp, 'T_read', T_read, 'S_max_required', S_max_required, ...
    'theta_cones', theta_cones, 'N_per_cone', N_per_cone, 'N_total', N_total, ...
    'gx_test', gx_test, 'gy_test', gy_test, 'gz_test', gz_test, 'wfm_test', wfm_test, ...
    'theta_test', theta_test, 'ktraj', ktraj, 't_ktraj', t_ktraj, ...
    'ktraj_adc', ktraj_adc, 't_adc', t_adc, 'slew_diag', slew_diag);

log{end+1} = sprintf('Total compile time: %.2f s', toc(tic0));
out.log = log;

end