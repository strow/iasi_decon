%
% NAME
%   kc2iasi - convolve kcarta to iasi channel radiances
%
% SYNOPSIS
%   [rad3, frq3] = kc2iasi(rkc, vkc)
%
% INPUTS
%   rkc     - kcarta grid radiances
%   vkc     - frequency grid for rkc
%
% OUTPUTS
%   rad3    - iasi radiances
%   frq3    - radiance frequency grid
%
% DISCUSSION
%   N    - kcarta to interferogram transform size
%   N3   - IASI user grid half-path size, opd / dx;
%   igm1 - N+1 pt single-sided high res igm from kcarta radiances
%   rad3 - IASI radiances from igm1
%   frq3 - IASI frequencies for rad3
% 
% HM, 17 Jul 2014
%

function [rad3, frq3] = kc2iasi(rkc, vkc)

%-----------------------------------
% set up interferometric parameters
%-----------------------------------

% kcarta params
dvk  = 0.0025;       % kcarta dv
Lmax = 1 / (2*dvk);  % kcarta nominal OPD

% IASI params
v1 = 645;            % iasi band low
v2 = 2760;           % iasi band high
vr = 20;             % out-of-band rolloff
dvc = 0.25;          % IASI dv
opd = 1 / (2*dvc);   % IASI OPD

% check that opd divides Lmax
% Lrat = Lmax / opd;
% if Lrat ~= floor(Lrat)
%   error('opd must divide Lmax')
% end

% find the transform size
% for k = 11 : 16
%   N = 2^k * Lrat;
%   if N * dvk >= v2, break, end
% end            

% set the transform size
N = 1638400;

% get Vmax and dx
Vmax = N * dvk;
dx = Lmax / N;
% isequal(Vmax, 1/(2*dx))

% single-sided undecimated interferogram size
N3 = opd / dx;
% isequal(N3, floor(N3))

%--------------------------------
% transform to channel radiances
%--------------------------------

% set kcarta radiance passband to the user grid
rkc = bandpass(vkc, rkc, v1, v2, vr);

% embed kcarta radiance in 0 to Vmax N+1 point grid
frq2 = (0:N)' * dvk;
rad2 = zeros(N+1, 1);
ix = interp1(frq2, (1:N+1)', vkc, 'nearest');
rad2(ix) = rkc;

% do the N+1 point cosine transform (as a 2*N point FFT)
igm1 = real(ifft([rad2; flipud(rad2(2:N,1))]));
igm1 = igm1(1:N+1,1);

% apply the IASI apodization
dtmp = (0:N3)' * dx;
apod = gaussapod(dtmp, 2);
igm1(1:N3+1) = igm1(1:N3+1) .* apod;

% truncate the single-sided high res igm1 to the instrument 
% OPD, and transform back to radiance
rad3 = real(fft([igm1(1:N3+1,1); flipud(igm1(2:N3,1))]));
frq3 = (0:N3)' * dvc;

% return just the IASI band
ix = interp1(frq3, (1:N3+1)', v1:dvc:v2, 'nearest');
rad3 = rad3(ix);
frq3 = frq3(ix);

