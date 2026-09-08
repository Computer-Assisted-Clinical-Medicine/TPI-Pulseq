% ==========================================================================
% TPI_compiler.m
% Headless (non-GUI) parameter script for the TPI Pulseq compiler.
%
% Use this instead of TPI_UI.m when running on a machine without a
% display, in a batch job, or when scripting parameter sweeps. Sets
% every field required by tpi_compile_advanced.m explicitly, then
% compiles and (optionally) saves the .seq file and diagnostic plots.
%
% Usage:
%   1. Edit the parameter block below to the desired protocol.
%   2. Run the script. out.log is printed to the console.
%   3. Adjust save_seq / save_plots / output_dir as needed.
%
% Dependencies:
%   - Pulseq (mr.*)                    github.com/pulseq/pulseq
%   - calc_cone_distribution.m, design_tpi_projection.m,
%     compute_slew_vs_cone_angle.m, tpi_compile_advanced.m (this repo)
%
% Author:  Valentin Jost
% Date:    2026-09
% ==========================================================================

clear; clc;

%% ---- Parameters (mirrors every field set by TPI_UI.m) ----
prm.scanner            = 'Skyra';     % 'Skyra' | 'Trio' | 'TerraW60'
prm.nuclei             = 'Na';        % 'Na' | 'H'

prm.fov                = 240e-3;      % field of view [m]
prm.Nx                 = 60;          % isotropic matrix size
prm.TR                 = 100e-3;      % repetition time [s]

prm.p                  = 0.4;         % twist onset fraction (k0 = p*kmax)
prm.Tread              = 7.6e-3;      % target effective readout time [s]
prm.flip_angle_deg     = 90;          % excitation flip angle [deg]

prm.ramp_mode          = 'linear';    % 'linear' | 'smooth'
prm.ramp_stretch       = 1.0;         % |G| ramp duration margin (>= 1)

prm.transition_mode    = 'none';      % 'none' | 'taper'
prm.transition_stretch = 1.0;         % taper duration multiplier

prm.avoid_equator      = false;       % exclude exact equatorial ring (odd NR only)
prm.symmetric_cones    = true;        % enforce mirrored ring pair symmetry

prm.ro_os              = 1.5;         % readout oversampling factor
prm.TEx                = 1;           % 1 or 2: which echo/contrast to compile
prm.TE2_delay          = 4.5e-3;      % [s], only used when TEx == 2
prm.undersample_factor = 1;           % additional undersampling divisor (>= 1)

%% ---- Output options ----
save_seq   = true;
save_plots = true;
output_dir = fullfile(pwd, 'output');
seq_name   = sprintf('TPI_%dmm_Nx%d_p%02d_%s_TE%d', ...
    round(prm.fov*1e3), prm.Nx, round(prm.p*10), prm.scanner, prm.TEx);

if (save_seq || save_plots) && ~isfolder(output_dir)
    mkdir(output_dir);
end

%% ---- Compile ----
fprintf('Compiling TPI sequence (headless) ...\n');
out = tpi_compile_advanced(prm);
fprintf('%s\n', strjoin(out.log, newline));

%% ---- Save .seq file ----
if save_seq
    seqpath = fullfile(output_dir, [seq_name '.seq']);
    out.seq.write(seqpath);
    fprintf('Sequence saved: %s\n', seqpath);
end

%% ---- Save diagnostic plots (k-space trajectory + slew sweep) ----
if save_plots
    fig1 = figure('Name', 'Full 3D k-space trajectory', 'Visible', 'off');
    plot3(out.ktraj(1,:), out.ktraj(2,:), out.ktraj(3,:), 'b-'); hold on;
    plot3(out.ktraj_adc(1,:), out.ktraj_adc(2,:), out.ktraj_adc(3,:), 'r.', 'MarkerSize', 2);
    axis equal; grid on; box on;
    xlabel('kx'); ylabel('ky'); zlabel('kz');
    title(sprintf('TPI k-space - %d projections, p=%.2f', out.N_total, prm.p));
    saveas(fig1, fullfile(output_dir, [seq_name '_ktraj.png']));
    close(fig1);

    fig2 = figure('Name', 'Slew rate vs. cone angle', 'Visible', 'off');
    d = out.slew_diag;
    plot(d.theta_deg, d.Snorm, 'm-', 'LineWidth', 1.6); hold on;
    plot(d.theta_deg, d.Saxis, 'k--', 'LineWidth', 1.2);
    plot(d.theta_deg, d.Sanalytical, 'c-.', 'LineWidth', 1.4);
    yline(d.limit, 'k:', sprintf('Hardware limit (%.0f T/m/s)', d.limit));
    xline(90, 'g:', 'kz = 0 (equator)');
    legend({'S_{vector}','S_{per-axis max}','S_{analytical} (Boada Eq.11)'}, 'Location', 'best');
    xlabel('Cone polar angle \theta_0 [deg]'); ylabel('Peak slew rate [T/m/s]');
    title('Slew rate vs. cone angle');
    grid on; box on;
    saveas(fig2, fullfile(output_dir, [seq_name '_slew.png']));
    close(fig2);

    fprintf('Diagnostic plots saved to: %s\n', output_dir);
end