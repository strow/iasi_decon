%
% NAME
%   iasi_decon - deconvolve IASI channel radiances
%
% SYNOPSIS
%   [rad2, frq2] = iasi_decon(rad1, frq1, dv2)
%
% INPUTS
%   rad1   - IASI channel radiances, m x n array
%   frq1   - IASI channel frequencies, m-vector
%   dv2    - optional deconvolution grid step, default 0.1 1/cm
%
% OUTPUTS
%   rad2   - deconvolved radiances, k x n array
%   frq2   - deconvolution frequency grid, k-vector
%
% DISCUSSION
%   see doc/finterp.pdf for the derivations used here.
%
%   iasi_decon calls isclose.m from airs_decon/test, bandpass.m 
%   from ccast/source, and gaussapod.m from /asl/matlib/fconv
%
%   iasi_decon and kc2iasi are similar, but kc2iasi convolves 
%   kcarta radiance to iasi channels while iasi_decon deconvolves
%   iasi channel radiances to an intermediate grid, typically at 
%   a 0.1 1/cm spacing.
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
  dv2 = 0.1;        % outut freq step
end

% IASI params
v1 = 645;            % iasi band low
v2 = 2760;           % iasi band high
vr = 20;             % out-of-band rolloff
vb = v2 + vr;        % transform max
dv1 = 0.25;          % IASI dv

% get rational approx to dv1/dv2
[m1, m2] = rat(dv1/dv2);
if ~isclose(m1/m2, dv1/dv2, 4)
  error('no rational approximation for dv1 / dv2')
end

% get the tranform sizes
for k = 4 : 24
  if m2 * 2^k * dv1 >= vb, break, end
end
N1 = m2 * 2^k;
N2 = m1 * 2^k;

% get (and check) dx
dx1 = 1 / (2*dv1*N1);
dx2 = 1 / (2*dv2*N2);
if ~isclose(dx1, dx2, 4)
  error('dx1 and dx2 are different')
end
dx = dx1;

% fprintf(1, 'iasi_decon: N1 = %7d, N2 = %5d, dx = %6.3e\n', N1, N2, dx);

%-------------------------------
% deconvolve the IASI radiances
%-------------------------------

% add a 5 cm-1 in-band rolloff at band edges
rad1 = bandpass(frq1, rad1, v1+5, v2-5, 5);

% embed IASI radiances in a 0 to Vmax grid
ftmp = (0:N1)' * dv1;
rtmp = zeros(N1+1, nobs);
[ix, jx] = seq_match(ftmp, frq1);
rtmp(ix, :) = rad1(jx, :);

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
ix = find(v1 <= frq2 & frq2 <= v2);
rad2 = rad2(ix, :);
frq2 = frq2(ix);

