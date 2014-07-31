%
% NAME
%   iasi_params - IASI interferometric parameters
%
% SYNOPSIS
%   user = iasi_params
%
% OUTPUT FIELDS
%   v1    - band low
%   v2    - band high
%   vr    - out-of band rolloff
%   dv    - frequency step, 1/(2*opd)
%   opd   - max optical path difference
%   npts  - point in frequency grid
%   freq  - IASI frequency grid
%
% DISCUSSION
%   field names follow the CrIS user-grid definitios
%
% AUTHOR
%   H. Motteler, 30 July 2014
%

function user = iasi_params

user.v1 = 645;
user.v2 = 2760;
user.dv = 0.25;
user.opd = 2.0;
user.vr = 20;

user.npts = round((user.v2 - user.v1) / user.dv) + 1;
user.freq = user.v1 + (0 : user.npts-1)' * user.dv;

