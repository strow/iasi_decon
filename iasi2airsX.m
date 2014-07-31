%
% NAME
%   iasi2airsX -- translate IASI to AIRS radiances
%
% SYNOPSIS
%   [arad, afrq] = iasi2airsX(rad1, frq1, sfile, cfreq, dvb)
%
% INPUTS
%   rad1   - IASI channel radiances, m x k array
%   frq1   - IASI channel frequencies, m-vector
%   sfile  - AIRS HDF4 SRF tabulation file
%   cfrq   - desired AIRS channel frequencies
%   dvb    - optional deconvolution grid step
%
% OUTPUTS
%   arad   - simulated AIRS radiances, n x k array
%   afrq   - AIRS channel frequencies, n-vector
%
% DISCUSSION
%   iasi2airsX translates IASI to AIRS radiances by interpolating
%   IASI to an intermediate grid and then convolving with the AIRS
%   channel response functions
%
%   iasi2airsX is used to test interpolation vs deconvolution to 
%   an intermediate grid followed by convolution to AIRS channels
%
% AUTHOR
%   H. Motteler, 20 July 2014
%

function [arad, afrq] = iasi2airsX(irad, ifrq, sfile, cfrq, dvb)

% check that array sizes match
ifrq = ifrq(:);
[m, nobs] = size(irad);
if m ~= length(ifrq)
  error('irad and ifrq sizes do not match')
end

% set default deconvolution grid step
if nargin < 5
  dvb = 0.05;
end

% IASI direct interpolation to AIRS
frq2 = (ifrq(1) : dvb : ifrq(end))';
rad2 = interp1(ifrq, irad, frq2, 'spline', 'extrap');

% get the AIRS convolution matrix
[sconv, sfrq, afrq] = mksconv1(sfile, cfrq, dvb);

% reconvolve to AIRS radiances
[ix, jx] = seq_match(sfrq, frq2);
rtmp = zeros(length(sfrq), nobs);
rtmp(ix, :) = rad2(jx, :);
arad = sconv * rtmp;

