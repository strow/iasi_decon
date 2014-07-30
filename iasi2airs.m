%
% NAME
%   iasi2airs -- translate IASI to AIRS radiances
%
% SYNOPSIS
%   [arad, afrq] = iasi2airs(rad1, frq1, sfile, cfreq, dvb)
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
%   iasi2airs translates IASI to AIRS radiances by deconvolving
%   IASI to an intermediate grid and then reconvolving with the 
%   AIRS channel response functions
%
%   note: the frequencies returned in afrq are from the AIRS SRF
%   tabulation.  These are matched in mksconv2 with the supplied
%   values (from cfreq) and returned if they are within 0.04 1/cm.
%
% COPYRIGHT
%   Copyright 2012-2014, Atmospheric Spectroscopy Laboratory.  
%   This code is distributed under the terms of the GNU GPL v3.
%
% AUTHOR
%   H. Motteler, 20 July 2014
%

function [arad, afrq] = iasi2airs(irad, ifrq, sfile, cfrq, dvb)

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

% deconvolve the IASI radiances
[rad2, frq2] = iasi_decon(irad, ifrq, dvb);

% get the AIRS convolution matrix
[sconv, sfrq, afrq] = mksconv2(sfile, cfrq, dvb);

% reconvolve to AIRS radiances
[ix, jx] = seq_match(sfrq, frq2);
rtmp = zeros(length(sfrq), nobs);
rtmp(ix, :) = rad2(jx, :);
arad = sconv * rtmp;

