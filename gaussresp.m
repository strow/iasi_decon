function r = gaussresp(v, L)

% function r = gaussresp(v, L)
%
% truncated-gaussian response function, via cosine transfrom
% of a truncated gaussian apodization
%
% inputs
%   v - wavenumbers (equally spaced, ascending from 0)
%   L - max path length
%
% output
%   r - response function of v

if nargin == 1
  L = 1;
end

L1 = L;
v = v(:);
n = length(v);

cospts = 2^nextpow2(n) + 1;
dv = v(2) - v(1);
vmax = dv * (cospts-1);
dd = 1/(2*vmax);
L1ind = round(L1/dd) + 1;
L1a = dd * (L1ind-1);
L1pts = 0:dd:L1a;

intf = zeros(cospts,1);
intf(1:L1ind) = gaussapod(L1pts', L1);
spec = real(fft([intf; flipud(intf(2:cospts-1,1))]));
spec = spec(1:cospts) / max(spec(1:cospts));

r = spec(1:n);

