function r = rms(A)

% function r = rms(A)
%
% returns RMS average of vector or matrix A

[m,n] = size(A);

r = norm(A(:)) / sqrt(m*n);

