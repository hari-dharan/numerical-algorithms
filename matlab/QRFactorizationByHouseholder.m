function [Q, R] = QRFactorizationByHouseholder(A)
%QRFACTORIZATIONBYHOUSEHOLDER Full QR Factorization using Householder transformations.
[m,n] = size(A);
Q = Householder(A(:,1), 0);
R = Q*A;
for i=2:n
    HH = Householder(R(i:m,i), i-1);
    Q = Q*HH;
    R = HH*R;
end
end