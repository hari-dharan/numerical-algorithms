function [Q, R] = QRFactorizationByGivens(A)
%QRFACTORIZATIONBYGIVENS Full QR Factorization using Givens rotations. 
[m,n] = size(A);
R = A;
Q = eye(m);
for j=1:n
    for i=j+1:m
        G = Givens(R(:,j), j, i);
        R = G*R;
        Q = Q*G';
    end
end
end