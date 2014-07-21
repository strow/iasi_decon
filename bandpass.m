%
% NAME
%   bandpass  -- bandpass filter
%
% SYNOPSIS
%   function dout = bandpass(vin, din, v1, v2, vr);
%
% INPUTS
%   vin   - input frequency grid, m-vector
%   din   - input data, m x n array, column order
%   v1    - passband start
%   v2    - passband end
%   vr    - optional rolloff width
%
% OUTPUT
%   dout  - din rolled off outside of [v1, v2]
%
% NOTES
%   the rolloff is a cosine, fit to the rolloff width
%
%   if the rolloff width vr is not specified, the rolloff is from
%   vin(1) to v1 and from v2 to vin(end).  If vr is specified and
%   v1 - vr < vin(1), then vin(1) is the low end of the rolloff,
%   similarly if vin(end) < v2 + vr then vin(end) is the high end.
%
%   the passband is taken as the smallest interval at the vin grid
%   points that spans [v1, v2].  
%
% AUTHOR
%    H. Motteler, 20 Apr 2012
%

function dout = bandpass(vin, din, v1, v2, vr)

% make vin a column vector
vin = vin(:);

% check that inputs vin and din conform
[m, n] = size(din);
if length(vin) ~= m
  error('vin length and din rows differ')
end

% set vr default value
if nargin < 5
  vr = vin(m);
end

% get passband indices in vin
j1 = max(find(vin <= v1));
j2 = min(find(v2 <= vin));

% get indices for filter wings 
k1 = max([find(vin <= v1 - vr); 1]);
k2 = min([find(v2 + vr <= vin); m]);

% get sizes of each segment
n1 = j1 - k1;      % points in LHS rolloff
n2 = k2 - j2;      % points in RHS rolloff
n3 = j2 - j1 + 1;  % points in passband

% scale cosine for the rolloffs
f1 = (1+cos(pi+(0:n1-1)*pi/n1))/2;
f2 = (1+cos((1:n2)*pi/n2))/2;

% build the filter
filt = zeros(m, 1);
filt(k1:k2, 1) = [f1, ones(1,n3), f2]';

% apply the filter
for i = 1 : n
  dout(:,i) = din(:,i) .* filt;
end

