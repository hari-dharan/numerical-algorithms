function [eigvec, eigval] = RayleighQuotientIteration(A, v, tol)
%RAYLEIGHQUOTIENTITERATION Find eigenvector and eigenvalue of A.
%   Rayleigh Quotient Iteration returns the smallest eigenvalue and corresponding eigenvector.
%   The iteration is defined by v_k1 = (A-sI)v_k, where s is a shift and I
%   is the identity. The algorithm instead solve (A-sI)v_k1 = v_k to avoid
%   computing the inverse. Since the shift is updated at every step, LU
%   Factorization is not used.
%   This iteration runs while the error between updates is greater than tol.
%   v_k is v at k-th iteration and v_k1 is v at k+1-th iteration.
v_k = zeros(size(v));
v_k1 = v;
while (norm(v_k1 - v_k) / norm(v_k)) > tol
    v_k = v_k1;
    sigma = (v_k1'*A*v_k1) / (v_k1'*v_k1);
    w = linsolve(A-sigma*eye(size(A)), v_k1);
    v_k1 = w/max(w);
end
sigma = (v_k1'*A*v_k1) / (v_k1'*v_k1);
eigvec = v_k1;
eigval = sigma;
end