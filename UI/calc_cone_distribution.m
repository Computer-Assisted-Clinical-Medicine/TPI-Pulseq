function [theta_cones, N_per_cone, N_total] = calc_cone_distribution(kmax, fov, p, undersample_factor)
% ==========================================================================
% CALC_CONE_DISTRIBUTION  Distribute TPI projections across polar cones.
%
% Implements the ring/cone distribution of Boada et al. (1997), Eq. [1].
% The k-space sphere is divided into NR rings (cones of constant polar
% angle theta), spaced by 1/FOV to satisfy the angular Nyquist criterion.
% Ring spacing depends only on angular coverage, not on the twist factor
% p (the twist changes azimuthal, not polar, sampling density).
%
% Projections per ring follow the conventional-PI density, reduced by:
%   p                   - TPI twist efficiency (Boada, p.709): a fixed
%                          gain tied to the trajectory design.
%   undersample_factor  - additional user-chosen undersampling below
%                          Nyquist (independent of trajectory physics).
%
% Inputs:
%   kmax                - maximum k-space radius [1/m]
%   fov                 - field of view [m]
%   p                   - twist onset fraction, 0 < p <= 1
%   undersample_factor  - undersampling divisor, >= 1 (default 1)
%
% Outputs:
%   theta_cones - polar angles per cone [rad], length NR
%   N_per_cone  - number of projections per cone
%   N_total     - total projection count
%
% Author:  Valentin Jost
% Date:    2026-09
% ==========================================================================

if nargin < 4 || isempty(undersample_factor)
    undersample_factor = 1;
end

% --- Ring (polar) sampling: angular Nyquist at kmax ---
NR = round(pi * kmax * fov / undersample_factor);
NR = max(NR, 2);
NR = NR + mod(NR, 2);              % force even NR (Boada Eq. [1])

n = (0:NR-1);
theta_cones = (n + 0.5) / NR * pi; % pole-safe ring centers

% --- Azimuthal (spoke) sampling per ring ---
N_per_cone_PI = 2*pi*kmax*fov .* sin(theta_cones);          % conventional PI baseline
N_per_cone    = round(p * N_per_cone_PI / undersample_factor);
N_per_cone    = max(N_per_cone, 1);                          % at least one spoke per ring

N_total = sum(N_per_cone);

end