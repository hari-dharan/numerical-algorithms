function [Q, R] = QRFactorization(A)
%QRFACTORIZATION Reduced QR Factorization using Modified Gram-Schmidt.
%   A = QR, where A is (nxm), Q is (nxm), R is (mxm) and n>m.
m = size(A,2);
R = zeros(m,m);
for i=1:m
    for j=1:i-1
        R(j,i) = A(:,i)'*A(:,j);
        A(:,i) = A(:,i) - R(j,i)*A(:,j);
    end
    R(i,i) = norm(A(:,i));
    A(:,i) = A(:,i)/R(i,i);
end
Q = A;
end