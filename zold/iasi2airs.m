%
% NAME
%   iasi2airs -- translate IASI to AIRS radiances
%
% SYNOPSIS
%   [arad, afrq] = iasi2airs(irad, ifrq, sfile, cfrq, dvb)
%
% INPUTS
%   irad   - IASI channel radiances, m x k array
%   ifrq   - IASI channel frequencies, m-vector
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
%   cfrq is desired channel frequencies and afrq is nominal center
%   frequencies of the tabulated SRFs.  If these differ by more than
%   0.02 1/cm, a warning is given.  iasi2airs does not sort afrq by
%   frequency.
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
  dvb = 0.1;
end

% deconvolve the IASI radiances
[rad2, frq2] = iasi_decon(irad, ifrq, dvb);

% get the AIRS convolution matrix
[sconv, sfrq, afrq] = mksconv1(sfile, cfrq, dvb);

% reconvolve to AIRS radiances
[ix, jx] = seq_match(sfrq, frq2);
rtmp = zeros(length(sfrq), nobs);
rtmp(ix, :) = rad2(jx, :);
arad = sconv * rtmp;

