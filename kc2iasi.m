%
% NAME
%   kc2iasi - convolve kcarta to IASI channel radiances
%
% SYNOPSIS
%   [rad2, frq2] = kc2iasi(rad1, frq1)
%
% INPUTS
%   rad1    - kcarta grid radiances
%   frq1    - frequency grid for rad1
%
% OUTPUTS
%   rad2    - IASI radiances
%   frq2    - IASI frequency grid
%
% DISCUSSION
%   see finterp.pdf
% 
% HM, 17 Jul 2014
%

function [rad2, frq2] = kc2iasi(rad1, frq1)

%-----------------------------------
% set up interferometric parameters
%-----------------------------------

% kcarta param
dv1  = 0.0025;       % kcarta dv

% IASI params
v1 = 645;            % iasi band low
v2 = 2760;           % iasi band high
vr = 20;             % out-of-band rolloff
dv2 = 0.25;          % IASI dv

% get rational approx to dv1/dv2
[m1, m2] = rat(dv1/dv2);
if ~isclose(m1/m2, dv1/dv2, 4)
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
% take kcarta to IASI radiances
%-------------------------------

% set kcarta radiance passband to the user grid
rad1 = bandpass(frq1, rad1, v1, v2, vr);

% embed kcarta radiance in a 0 to Vmax grid
ftmp = (0:N1)' * dv1;
rtmp = zeros(N1+1, 1);
ix = interp1(ftmp, (1:N1+1)', frq1, 'nearest');
rtmp(ix) = rad1;

% radiance to interferogram
igm1 = real(ifft([rtmp; flipud(rtmp(2:N1,1))]));
igm1 = igm1(1:N1+1,1);

% apply the IASI apodization
dtmp = (0:N2)' * dx;
apod = gaussapod(dtmp, 2);
igm1(1:N2+1) = igm1(1:N2+1) .* apod;

% interferogram to radiance
rad2 = real(fft([igm1(1:N2+1,1); flipud(igm1(2:N2,1))]));
frq2 = (0:N2)' * dv2;

% return just the IASI band
ix = interp1(frq2, (1:N2+1)', v1:dv2:v2, 'nearest');
rad2 = rad2(ix);
frq2 = frq2(ix);

