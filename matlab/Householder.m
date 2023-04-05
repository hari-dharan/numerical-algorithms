function HH = Householder(x,s)
%HOUSEHOLDER Householder reflection transformation given a vector x.
%   The parameter s indicates the size of the padding on the top left with
%   the identity matrix. When s = 0, the usual Householder matrix for a (nx1)
%   vector x is returned with size (nxn). When s > 0, the Householder
%   matrix returned is of size ((n+s)x(n+s)).
L = length(x);
N = sqrt(x'*x);
x(1) = x(1) - N;
H = eye(L) - 2*((x*x')/(x'*x));
if s == 0
    HH = H;
else
    HH = [eye(s), zeros(s,L); zeros(L,s), H];
end