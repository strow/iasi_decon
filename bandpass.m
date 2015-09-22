%
% NAME
%   bandpass  -- bandpass filter with raised cosine wings
%
% SYNOPSIS
%   function b = bandpass(v, d, v1, v2, r1, r2);
%
% INPUTS
%   v   - frequency grid, m-vector
%   d   - input data, m x n array
%   v1  - passband start frequency
%   v2  - passband end frequency
%   r1  - optional LHS rolloff
%   r2  - optional RHS rolloff 
%
% OUTPUT
%   b   - d filtered, m x n array
%
% NOTES
%   if r1 is not specified the value is set to 10
%   if r2 is not specifed the value for r1 is used
%
%   this function differs from earlier versions in that it allows
%   for the separate specification of the LHS and RHS wings, and in
%   that the passband and wings are purely functions of frequency
%   rather than being adjusted to fall on grid points spanning the
%   passband.
%
% AUTHOR
%    H. Motteler, 15 Jun 2015
%

function b = bandpass(v, d, v1, v2, r1, r2)

% make v a column vector
v = v(:);

% check that v and d conform
[m, n] = size(d);
if length(v) ~= m
  error('v and d do not conform')
end

% set rolloff defaults
if nargin == 4
  r1 = 10;
  r2 = 10;
elseif nargin == 5
  r2 = r1;
end

% quietly fix long r1 or r2
vL = max(v(1), v1 - r1);
vH = min(v2 + r2, v(end));

% partition indices
c2 = vL <= v & v < v1;
c3 = v1 <= v & v < v2;
c4 = v2 <= v & v < vH;

% build the filter
f = zeros(m, 1);
f(c2) = (1 + cos(-pi * (v1 - v(c2)) / (v1 - vL))) / 2;
f(c3) = 1;
f(c4) = (1 + cos( pi * (v2 - v(c4)) / (v2 - vH))) / 2;

% apply the filter
for i = 1 : n
  b(:,i) = d(:,i) .* f;
end

