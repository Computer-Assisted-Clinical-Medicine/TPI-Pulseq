function [SNRmap, SNR0map, sigma] = computeSNRmap(vol4D, mask, noiseType)
if nargin < 3, noiseType = 'gaussian'; end

spatial_shape = size(mask);              % always true [Nx Ny Nz]
nvox_spatial  = numel(mask);
NTE           = numel(vol4D) / nvox_spatial;   % robust even if singleton dim was dropped

if mod(NTE,1) ~= 0
    error('computeSNRmap: numel(vol4D) not divisible by numel(mask) — check inputs.');
end
NTE = round(NTE);

data_flat = reshape(vol4D, nvox_spatial, NTE);  % force correct shape explicitly
mask_flat = reshape(mask, nvox_spatial, 1);

switch lower(noiseType)
    case 'gaussian'
        sigma = estimateGaussianNoise(vol4D, mask);
    case 'rician'
        sigma = estimateRicianNoise(vol4D, mask);
    otherwise
        error('computeSNRmap: noiseType must be ''gaussian'' or ''rician''.');
end

SNR_flat = NaN(nvox_spatial, NTE);

for iTE = 1:NTE
    M = data_flat(:, iTE);
    valid = mask_flat & isfinite(M) & M > 0;

    switch lower(noiseType)
        case 'gaussian'
            SNR_flat(valid, iTE) = M(valid) / sigma(iTE);
        case 'rician'
            A_hat = sqrt(max(M(valid).^2 - sigma(iTE)^2, 0));
            SNR_flat(valid, iTE) = A_hat / sigma(iTE);
    end
end

SNR0_flat = SNR_flat(:, 1);

if NTE > 1
    SNRmap = reshape(SNR_flat, [spatial_shape, NTE]);
else
    SNRmap = reshape(SNR_flat, spatial_shape);   % avoid trailing singleton mismatch on output too
end
SNR0map = reshape(SNR0_flat, spatial_shape);

fprintf('computeSNRmap (%s): sigma range %.4e to %.4e | median SNR0 = %.2f\n', ...
    noiseType, min(sigma), max(sigma), median(SNR0_flat(isfinite(SNR0_flat))));
end