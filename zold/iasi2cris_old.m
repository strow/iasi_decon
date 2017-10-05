%
% NAME
%   iasi2cris -- translate IASI to CrIS radiances
%
% SYNOPSIS
%   [crad, cfreq] = iasi2cris(rad1, frq1, opt1)
%
% INPUTS
%   rad1   - IASI channel radiances, m x n array
%   frq1   - IASI channel frequencies, m-vector
%   opt1   - optional input parameters
%
% opt1 fields
%   hapod    - 0 = no apod (default), 1 = Hamming
%   resmode  - 'lowres' (default), 'hires1', or 'hires2'
%   nguard   - number of guard channels, default is 0
%
% OUTPUTS
%   crad   - simulated CrIS radiances, n x k array
%   cfrq   - CrIS channel frequencies, n-vector
%
% DISCUSSION
%   iasi2cris translates IASI to CrIS radiances with three direct
%   convolutions, one for each CrIS band with the relevant bandpass
%   filtering applied to the IASI radiances before the transform.
%
% COPYRIGHT
%   Copyright 2013-2014, Atmospheric Spectroscopy Laboratory.  
%   This code is distributed under the terms of the GNU GPL v3.
%
% AUTHOR
%   H. Motteler, 27 Oct 2014
%

function [crad, cfrq] = iasi2cris(rad1, frq1, opt1)

% defaults
hapod = 0;      % no Hamming apodization
nguard = 0;     % no guard channels

% process input options
if nargin == 3
  if isfield(opt1, 'hapod'), hapod = opt1.hapod; end
  if isfield(opt1, 'nguard'), nguard = opt1.nguard; end
end

% IASI params
iasi.v1 = 645;            % iasi band low
iasi.v2 = 2760;           % iasi band high
iasi.dv = 0.25;           % IASI dv

% CrIS params
wlaser = 773.1301;        % nominal value
bstr{1} = 'LW';           % band by number
bstr{2} = 'MW'; 
bstr{3} = 'SW';

% initialize outputs
crad = []; cfrq = []; 
bpts = []; bv1 = []; bv2 = [];

% loop on CrIS bands
for bi = 1 : 3

  % get the CrIS user grid
  band = bstr{bi};
  [inst, user] = inst_params(band, wlaser, opt1);

  % get the passband plus usable filter wings
  tv1 = max(iasi.v1, user.v1 - user.vr);
  tv2 = min(iasi.v2, user.v2 + user.vr);
  tvr = min(user.v1 - tv1, tv2 - user.v2);

  % trim IASI data to the filter span 
  ix = find(tv1 <= frq1 & frq1 <= tv2);
  ftmp = frq1(ix); 
  rtmp = rad1(ix, :); 

  % apply the bandpass filter
  rtmp = bandpass(ftmp, rtmp, user.v1, user.v2, tvr);

  % convolve to the CrIS user grid
  [rtmp, ftmp] = iasi_decon(rtmp, ftmp, user.dv, opt1);
  ftmp = ftmp(:);

  % option for hamming apodization
  if hapod
    rtmp = hamm_app(rtmp);
  end

  % return the CrIS user grid plus any guard channels
  dg = nguard * user.dv;
  ix = find(user.v1 - dg <= ftmp & ftmp <= user.v2 + dg);
  rtmp = rtmp(ix, :);
  ftmp = ftmp(ix);

  % concatenate output columns
  crad = [crad; rtmp];
  cfrq = [cfrq; ftmp];

end

