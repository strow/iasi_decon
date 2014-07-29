%
% NAME
%   iasi_decon - deconvolve iasi channel radiances
%
% SYNOPSIS
%   [rad2, frq2] = iasi_decon(rad1, frq1, dv2)
%
% INPUTS
%   rad1   - IASI channel radiances, m x n array
%   frq1   - IASI channel frequencies, m-vector
%   dv2    - optional deconvolution frequency step
%
% OUTPUTS
%   rad2   - deconvolved radiances, k x n array
%   frq2   - deconvolution frequency grid, k-vector
%
% DISCUSSION
%   see finterp.pdf
%
% HM, 17 Jul 2014
%

function [rad2, frq2] = iasi_decon(rad1, frq1, dv2)

% check that array sizes match
frq1 = frq1(:);
[m, nobs] = size(rad1);
if m ~= length(frq1)
  error('rad1 and frq1 sizes do not match')
end

%-----------------------------------
% set up interferometric parameters
%-----------------------------------

% default output frequency step
if nargin == 2
  dv2 = 0.05;        % outut freq step
end

% IASI params
v1 = 645;            % iasi band low
v2 = 2760;           % iasi band high
vr = 20;             % out-of-band rolloff
dv1 = 0.25;          % IASI dv

% get rational approx to dv1/dv2
[m1, m2] = rat(dv1/dv2);
if ~isclose(m1/m2, dv1/dv2)
  error('no rational approximation for dv1 / dv2')
end

% get the tranform sizes
for k = 5 : 24
  if m2 * 2^k * dv1 >= v2, break, end
end
N1 = m2 * 2^k;
N2 = m1 * 2^k;

% get (and check) dx
dx1 = 1 / (2*dv1*N1);
dx2 = 1 / (2*dv2*N2);
if ~isclose(dx1, dx2)
  error('dx1 and dx2 are different')
end
dx = dx1;

%-------------------------------
% deconvolve the IASI radiances
%-------------------------------

% embed IASI radiances in a 0 to Vmax grid
ftmp = (0:N1)' * dv1;
rtmp = zeros(N1+1, nobs);
ix = interp1(ftmp, (1:N1+1)', frq1, 'nearest');
rtmp(ix, :) = rad1;

% radiance to interferogram
igm1 = real(ifft([rtmp; flipud(rtmp(2:N1, :))]));
igm1 = igm1(1:N1+1, :);

% apply the inverse IASI apodization
dtmp = (0:N1)' * dx;
apod = gaussapod(dtmp, 2) * ones(1, nobs);
igm1 = igm1 ./ apod;

% extend the apodized interferogram
igm2 = zeros(N2+1, nobs);
igm2(1:N1+1, :) = igm1;

% interferogram to radiance
rad2 = real(fft([igm2(1:N2+1,:); flipud(igm2(2:N2,:))]));
frq2 = (0:N2)' * dv2;

% return just the IASI band
ix = interp1(frq2, (1:N2+1)', v1:dv2:v2, 'nearest');
rad2 = rad2(ix, :);
frq2 = frq2(ix);

