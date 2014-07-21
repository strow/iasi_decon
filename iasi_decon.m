%
% NAME
%   iasi_decon - deconvolve iasi channel radiances
%
% SYNOPSIS
%   [brad, bfrq] = iasi_decon(irad, ifrq, dvb)
%
% INPUTS
%   irad   - IASI channel radiances, m x k array
%   ifrq   - IASI channel frequencies, m-vector
%   dvb    - optional deconvolution frequency step
%
% OUTPUTS
%   brad   - deconvolved radiances
%   bfrq   - deconvolution frequency grid
%
% DISCUSSION
%   N    - kcarta to interferogram transform size
%   N3   - IASI user grid half-path size, opd / dx;
%   igm1 - N+1 pt single-sided high res igm from kcarta radiances
%
% HM, 17 Jul 2014
%

function [brad, bfrq] = iasi_decon(irad, ifrq, dvb)

%-----------------------------------
% set up interferometric parameters
%-----------------------------------

% default deconvolution grid
if nargin == 2
  dvb = 0.05;        % dcon freq step
end
Lmax = 1 / (2*dvb);  % transform OPD

% IASI params
v1 = 645;            % iasi band low
v2 = 2760;           % iasi band high
vr = 20;             % out-of-band rolloff
dvc = 0.25;          % IASI dv
opd = 1 / (2*dvc);   % IASI OPD

% check that opd divides Lmax
Lrat = Lmax / opd;
if Lrat ~= floor(Lrat)
  error('opd must divide Lmax')
end

% find the transform size
for k = 4 : 20
  N = 2^k * Lrat;
  if N * dvb >= v2, break, end
end            

% get Vmax and dx
Vmax = N * dvb;
dx = Lmax / N;
isequal(Vmax, 1/(2*dx))

% single-sided undecimated interferogram size
N3 = opd / dx;
isequal(N3, floor(N3))

%--------------------------------
% transform to channel radiances
%--------------------------------

% embed IASI radiances in a 0 to Vmax N3+1 point grid
frq1 = (0:N3)' * dvc;
rad1 = zeros(N3+1, 1);
ix = interp1(frq1, (1:N3+1)', ifrq, 'nearest');
rad1(ix) = irad;

% do an N3+1 point cosine transform (as a 2*N3 point FFT)
igm1 = real(ifft([rad1; flipud(rad1(2:N3,1))]));
igm1 = igm1(1:N3+1,1);

% apply the inverse IASI apodization
dtmp = (0:N3)' * dx;
apod = gaussapod(dtmp, 2);
igm1 = igm1 ./ apod;

% extend the apodized interferogram
igm2 = zeros(N+1, 1);
igm2(1:N3+1,1) = igm1;

% transform back to radiance
rad2 = real(fft([igm2(1:N+1,1); flipud(igm2(2:N,1))]));
frq2 = (0:N)' * dvb;

% return just the IASI band
ix = interp1(frq2, (1:N+1)', v1:dvb:v2, 'nearest');
brad = rad2(ix);
bfrq = frq2(ix);

