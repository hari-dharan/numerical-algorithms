function [eigvecs, eigvals] = QRIteration(A, k)
%QRITERATION QR Iteration for a real symmetric matrix A.
Q = eye(size(A));
for i=1:k
    [Qi, Ri] = qr(A);
    A = Ri*Qi;
    Q = Q*Qi;
end
eigvecs = Q;
eigvals = diag(A);
end