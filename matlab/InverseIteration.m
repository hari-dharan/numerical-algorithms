function [eigvec, eigval] = InverseIteration(A, v, sigma, tol)
%INVERSEITERATION Find eigenvector and eigenvalue of A.
%   Inverse Iteration returns the smallest eigenvalue and corresponding eigenvector.
%   The iteration is defined by v_k1 = (A-sI)v_k, where s is a shift and I
%   is the identity. The algorithm instead solve (A-sI)v_k1 = v_k to avoid
%   computing the inverse. This is solved using LU Factorization so that L
%   and U can be reused in the loop.
%   This iteration runs while the error between updates is greater than tol.
%   v_k is v at k-th iteration and v_k1 is v at k+1-th iteration.
[L, U] = lu(A-sigma*eye(size(A)));
v_k = zeros(size(v));
v_k1 = v;
while (norm(v_k1 - v_k) / norm(v_k)) > tol
    v_k = v_k1;
    x = linsolve(L, v_k1);
    w = linsolve(U, x);
    v_k1 = w/max(w);
end
eigvec = v_k1;
eigval = 1/max(w) + sigma;
end