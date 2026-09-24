% ==========================================================================
% TPI_Recon.m
% Reconstruction script for Twisted Projection Imaging (TPI) sodium MRI
%
% Sequence: single-readout-per-projection (no multi-echo), center-out TPI
%           trajectories as produced by tpi_compile_advanced.m / TPI_UI.m
%
% Pipeline:
%   1. Load .seq file, recompute trajectory, extract sequence definitions
%   2. Load raw Siemens .dat file via mapVBVD
%   3. Coil combination (if multi-channel array used)
%   4. Trajectory-vs-ADC sanity check (visual)
%   5. Reshape, normalize, and (optionally) average k-space data
%   6. NUFFT (density-compensated gridding / iNUFT) reconstruction
%   7. Display: regridded k-space, reconstructed image

% Dependencies:
%   - Pulseq MATLAB toolbox   (mr.Sequence)              github.com/pulseq/pulseq
%   - mapVBVD                 (Siemens raw data reader)
%   - nufft_3d                (non-uniform FFT)           github.com/marcsous/nufft_3d
%   - fft3c_new, ifft3c_new  (local utility functions)
%
% Author:  Valentin Jost
% Date:    2026-04
% ==========================================================================

clear; clc;

%% ---- Section 1: File Paths ----
% The .seq file and matching Siemens .dat file must share the same base name.
repo_root = fileparts(mfilename('fullpath'));
addpath(fullfile(repo_root, 'UI'));

% data_dir = fullfile(pwd);
data_dir = fullfile(repo_root,'ExampleData');
seqname  = 'TPI_240fov_2_5mm_p045_us3_TE1_smooth33';   % base name (no extension)

seqpath  = fullfile(data_dir, [seqname '.seq']);
datapath = fullfile(data_dir, [seqname '.dat']);

%% ---- Section 2: Load .seq File & Compute k-Space Trajectory ----
% seq.calculateKspacePP() integrates the gradient waveforms to give:
%   ktraj_adc  [3 x N_adc_total]  - k-space coordinates at every ADC sample
%   t_adc      [1 x N_adc_total]  - time stamps of ADC samples [s]
%   ktraj      [3 x N_grad_total] - full gradient-raster trajectory
%   t_ktraj    [1 x N_grad_total] - time stamps on gradient raster [s]
% N_adc_total = N_adc (samples/projection) x N_total (projections)

prevWarnState = warning('off', 'mr:restoreShape');

fprintf('Loading sequence file ... ');
seq = mr.Sequence();
seq.read(seqpath);
fprintf('done\n');

fprintf('Computing k-space trajectory ... ');
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation] = seq.calculateKspacePP();
fprintf('done\n');

warning(prevWarnState);   % restore whatever the warning state was before

kmax_meas = max(abs(ktraj_adc), [], 'all');   % [1/m], sanity check vs. stored kmax

%% ---- Section 3: Extract Sequence Parameters from .seq Definitions ----
% Requires tpi_compile_advanced.m to persist the following via
% seq.setDefinition(): 'FOV','Name','p','NR','N_total','kmax','G_amp',
% 'T_read','Nx','N_adc','ro_os'. No fallbacks are applied here - if a
% definition is missing, fix the compiler rather than patching the recon.
fprintf('Reading sequence definitions ... ');

required_defs = {'FOV','p','kmax','G_amp','T_read','N_total','NR','Nx','N_adc','ro_os'};
for i = 1:numel(required_defs)
    if isempty(seq.getDefinition(required_defs{i}))
        error(['TPI_Recon: required definition ''%s'' is missing from the .seq file.\n' ...
               'Update tpi_compile_advanced.m to call seq.setDefinition(''%s'', ...) before seq.write().'], ...
              required_defs{i}, required_defs{i});
    end
end

fov_vec = seq.getDefinition('FOV');
fov     = fov_vec(1);                       % isotropic FOV -> take first element
p       = seq.getDefinition('p');           % twist onset fraction (k0 = p*kmax)
kmax    = seq.getDefinition('kmax');        % max k-space radius [1/m]
G_amp   = seq.getDefinition('G_amp');       % gradient amplitude [T/m]
T_read  = seq.getDefinition('T_read');      % readout duration [s]
N_total = seq.getDefinition('N_total');     % total number of projections
NR      = seq.getDefinition('NR');          % number of cones
Nx      = seq.getDefinition('Nx');          % matrix size
N_adc   = seq.getDefinition('N_adc');       % samples per projection
ro_os   = seq.getDefinition('ro_os');       % readout oversampling factor

kspace_size = Nx;    % isotropic reconstruction grid size
navg        = 1;     % placeholder, re-derived from data in Section 7

fprintf('done\n');
fprintf('  FOV       = %.0f mm\n',    fov*1e3);
fprintf('  Nx        = %d\n',         Nx);
fprintf('  kmax      = %.1f 1/m\n',   kmax);
fprintf('  kmax_meas = %.1f 1/m\n',   kmax_meas);
fprintf('  p         = %.2f\n',       p);
fprintf('  N_total   = %d projections\n', N_total);
fprintf('  NR        = %d cones\n',   NR);
fprintf('  N_adc     = %d samples/projection\n', N_adc);
fprintf('  ro_os     = %.2f\n',       ro_os);
fprintf('  T_read    = %.3f ms\n',    T_read*1e3);

if abs(kmax_meas - kmax) / kmax > 0.01
    warning('kmax from trajectory (%.1f) deviates >1%% from stored kmax (%.1f) [1/m]', ...
        kmax_meas, kmax);
end
%% ---- Section 4: Load Raw Siemens Data (mapVBVD) ----
% mapVBVD returns a twix object; multi-RAID files return a cell array
% whose last element is the imaging scan.
%
% data_unsorted dims (mapVBVD .image.unsorted()):
%   dim 1: ADC samples  -> N_adc_total = N_adc * N_total
%   dim 2: coil channels
%   dim 3: repetitions/lines, Siemens acquisition order (NOT sorted)
% For TPI (single readout/projection, no phase encoding):
%   dim 3 = N_total                (no averages)
%         = N_total * navg         (multiple averages acquired)

fprintf('Loading raw data (mapVBVD) ... ');
twix_obj = mapVBVD(datapath);

if iscell(twix_obj)
    data_unsorted = twix_obj{end}.image.unsorted();   % last entry = imaging dataset
else
    data_unsorted = twix_obj.image.unsorted();
end
fprintf('done\n');
fprintf('  Raw data size: [%s]\n', num2str(size(data_unsorted)));   % expect [N_adc_total, n_coils, N_total*navg]

%% ---- Section 5: Coil Combination ----

n_coils = size(data_unsorted, 2);
fprintf('Detected coils: %d\n', n_coils);

if n_coils > 1
    fprintf('Coil combination ... ');
    data_unsorted = my_coilcombine2(data_unsorted, 2);

    fprintf('done\n');
end
% After combination: data_unsorted is [N_adc_total, 1, N_total*navg]

%% ---- Section 6: Trajectory Sanity Check (Visual) ----
% Overlay gradient trajectory, ADC sampling positions, and one raw FID
% to confirm timing alignment before gridding. TPI has a single readout
% window per projection (no multi-echo splitting needed here).

adc_len_total = size(data_unsorted, 1);
proj_plot = min(1000, N_total);   % representative projection, bounds-safe

figure('Name', 'TPI Trajectory Sanity Check');
hold on;
plot(t_ktraj, ktraj(1, :), 'r', 'DisplayName', 'kx');
plot(t_ktraj, ktraj(2, :), 'g', 'DisplayName', 'ky');
plot(t_ktraj, ktraj(3, :), 'b', 'DisplayName', 'kz');
plot(t_adc, ktraj_adc(3, :), 'ko', 'MarkerSize', 2, 'MarkerFaceColor', 'r', ...
    'DisplayName', 'kz ADC');
if adc_len_total >= proj_plot * N_adc
    fid_idx = (proj_plot-1)*N_adc + 1 : proj_plot*N_adc;
    plot(t_adc(fid_idx), abs(data_unsorted(fid_idx, 1, proj_plot)) * 1e5, ...
        'k--', 'DisplayName', sprintf('|FID| proj %d (x1e5)', proj_plot));
end
xlabel('Time [s]'); ylabel('k [1/m] / arb.');
title('TPI k-space trajectory vs. ADC windows');
legend('Location', 'best'); grid on;

TR_est = t_ktraj(2) - 0; % or read TR from definitions/estimate from RF spacing
tr_win = t_ktraj <= (min(t_excitation(2), t_ktraj(end)));  % first TR only
adc_win = t_adc <= min(t_excitation(2), t_adc(end));

figure('Name', 'TPI Trajectory Sanity Check (single TR)');
hold on;
plot(t_ktraj(tr_win), ktraj(1, tr_win), 'r', 'DisplayName', 'kx');
plot(t_ktraj(tr_win), ktraj(2, tr_win), 'g', 'DisplayName', 'ky');
plot(t_ktraj(tr_win), ktraj(3, tr_win), 'b', 'DisplayName', 'kz');
plot(t_adc(adc_win), ktraj_adc(3, adc_win), 'ko', 'MarkerSize', 4, ...
    'MarkerFaceColor', 'y', 'DisplayName', 'kz ADC samples');
xlabel('Time [s]'); ylabel('k [1/m]');
title('TPI trajectory vs. ADC windows — single TR');
legend('Location', 'best'); grid on;

%% ---- Section 7: Data Handling ----

% 7.1 Number of averages, inferred from data (not stored in .seq)
n_readouts_total = size(data_unsorted, 3);   % N_total * navg
assert(mod(n_readouts_total, N_total) == 0, ...
    ['Number of readouts (%d) is not divisible by N_total (%d).\n' ...
     'Scan may have been aborted or dummy scans are present.'], ...
    n_readouts_total, N_total);

navg = n_readouts_total / N_total;
fprintf('  N_adc            = %d  (samples/projection, from seq)\n', N_adc);
fprintf('  N_total          = %d  (projections,        from seq)\n', N_total);
fprintf('  n_readouts_total = %d  (from data)\n', n_readouts_total);
fprintf('  navg             = %d  (inferred: n_readouts / N_total)\n', navg);

% 7.2 Reshape raw data -> [N_adc, N_total, navg]
data = reshape(data_unsorted, N_adc, N_total, navg);
fprintf('  data shape       = [%s]  [N_adc, N_total, navg]\n', num2str(size(data)));

% 7.3 Reshape k-space trajectory -> [3, N_adc, N_total]
% (Identical for every average, so it is not replicated.)
kcoord = reshape(ktraj_adc, 3, N_adc, N_total);
fprintf('  kcoord shape     = [%s]  [3, N_adc, N_total]\n', num2str(size(kcoord)));

% ---- 7.4 Normalize k-space to nufft_3d's convention (cycles/FOV, Nyquist = ±Nx/2) ----
% nufft_3d.m's own docstring: "trajectory units are phase cycles/fov,
% Nyquist distance is 1 unit" -> om = k [1/m] * fov [m] = k/(2*kmax)*Nx
% Uses the stored kmax (not kmax_meas), since the measured peak slightly
% overshoots due to gradient-raster discretization.
kcoord_norm = kcoord / (2 * kmax) * Nx;   % = kcoord * fov, cycles/FOV convention

k_norm_max = max(abs(kcoord_norm(:)));
k_norm_limit = Nx / 2;  
fprintf('  |k_norm| max     = %.4f  (should be <= %.1f = Nx/2)\n', k_norm_max, k_norm_limit);
if k_norm_max > 1.001 * k_norm_limit
    warning('Normalized k-space exceeds Nx/2 (%.1f) - check kmax/Nx definitions in compiler.', k_norm_limit);
end

% 7.5 Optional k-space averaging (before NUFFT, for SNR)
k_avg = true;   % false keeps averages separate (e.g. for motion inspection)

if k_avg && navg > 1
    fprintf('  Averaging %d repetitions in k-space ... ', navg);
    data       = mean(data, 3);   % [N_adc, N_total, 1]
    navg_recon = 1;
    fprintf('done\n');
else
    navg_recon = navg;
end
data = squeeze(data);   % [N_adc, N_total] or [N_adc, N_total, navg]
fprintf('  Data shape after averaging: [%s]\n', num2str(size(data)));

% 7.6 Flatten data and coordinates for NUFFT input
%   kcoord_flat [3, N_adc*N_total]  normalized coordinates
%   data_flat   [N_adc*N_total, 1]  complex k-space values
kcoord_flat = reshape(kcoord_norm, 3, []);
data_flat   = reshape(data, [], 1);
fprintf('  NUFFT input: kcoord [%s], data [%s]\n', ...
    num2str(size(kcoord_flat)), num2str(size(data_flat)));

% 7.7 Diagnostic: mean radial k-space profile (T2* decay shape, twist
% onset at k0 = p*kmax, DC offset check)
k_radius     = sqrt(sum(kcoord(:,:,1).^2, 1));
mean_profile = mean(abs(reshape(data, N_adc, [])), 2);

[k_sorted, sort_idx] = sort(k_radius);
profile_sorted = mean_profile(sort_idx);

figure('Name', 'Mean radial k-space profile');
plot(k_sorted, profile_sorted, 'b');
xlabel('k-space radius [1/m]'); ylabel('Mean |signal| [a.u.]');
title('TPI radial signal profile (mean over all projections)');
xline(kmax*p, '--r', sprintf('k_0 = p*k_{max} = %.1f 1/m', kmax*p));
grid on;

%% ---- Section 8: NUFFT Reconstruction ----
% nufft_3d options: 'gpu' GPU accel, 'u' grid oversampling, 'J' kernel
% width [cells], 'radial' built-in radial density compensation.

% 8.1 Pre-allocate outputs
img         = zeros(kspace_size, kspace_size, kspace_size, navg_recon);
kspace_grid = zeros(kspace_size, kspace_size, kspace_size, navg_recon);

% 8.2 Build the NUFFT operator once (trajectory is identical per average)
fprintf('Building NUFFT operator ... ');
nufft_obj = nufft_3d(kcoord_flat, kspace_size, 'gpu', 0, 'u', 4, 'J', 6, 'radial', 0);
fprintf('done\n');

maxit  = 1;   % 1 = density-compensated gridding; >1 = conjugate-gradient iterations
damp   = 0;   % L2 regularization on ||x|| (0 = none, for gridding)
weight = 1;   % uniform sample weighting (replace with a DCF vector if available)

% 8.3 Reconstruction loop over averages
fprintf('Reconstruction ');
for avg_i = 1:navg_recon
    fprintf('(avg %d/%d) ... ', avg_i, navg_recon);

    if navg_recon == 1
        data_i = data_flat;
    else
        data_i = reshape(data(:, :, avg_i), [], 1);
    end

    im = nufft_obj.iNUFT(data_i, maxit, damp, weight);
    im = reshape(im, kspace_size, kspace_size, kspace_size);
    img(:,:,:,avg_i) = im;

    kspace_i = fft3c_new(im);                              % regrid for inspection
    kspace_grid(:,:,:,avg_i) = kspace_i;

    fprintf('done\n');
end

% Combine averages by magnitude (incoherent averaging, post-recon)
if navg_recon > 1
    img_avg         = mean(abs(img),         4);
else
    img_avg         = abs(img(:,:,:,1));
end

fprintf('Reconstruction complete.\n');

%% ---- Section 9: Display ----
N = kspace_size;

% 9.1 Regridded k-space, central slice (log magnitude)
% kspace_slice = squeeze(kspace_grid(:, N/2+1, :, 1));
kspace_slice = squeeze(kspace_grid(N/2+1, :, :, 1));
% kspace_slice = squeeze(kspace_grid(:, :, N/2+1, 1));

figure('Name', 'Regridded k-space (central slice)');
imagesc(log10(abs(kspace_slice)));
% clim([-12 -5]);
cb=colorbar; cb.Label.String = 'log(signal) [a.u.]'; axis image;
title(sprintf('TPI regridded k-space - %d projections, p=%.2f', N_total, p));
xlabel('kx'); ylabel('kz');

% 9.2 Reconstructed image, three orthogonal centre slices
figure('Name', 'TPI reconstructed image');

subplot(1,3,1); imagesc(abs(img_avg(:,:,N/2+1)));
axis image off; colormap gray; colorbar; title('Axial (z = centre)');

subplot(1,3,2); imagesc(squeeze(abs(img_avg(:,N/2+1,:))));
axis image off; colorbar; title('Coronal (y = centre)');

subplot(1,3,3); imagesc(squeeze(abs(img_avg(N/2+1,:,:))));
axis image off; colorbar; title('Sagittal (x = centre)');

sgtitle(sprintf('TPI Na - FOV=%.0fmm, Nx=%d, p=%.2f', fov*1e3, Nx, p));