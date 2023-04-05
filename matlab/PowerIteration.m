function [eigvec, eigval] = PowerIteration(A, v, tol)
%POWERITERATION Find eigenvector and eigenvalue of A.
%   Power Iteration returns the largest eigenvalue and corresponding eigenvector.
%   This iteration runs while the error between updates is greater than tol.
%   v_k is v at k-th iteration and v_k1 is v at k+1-th iteration.
v_k = zeros(size(v));
v_k1 = v;
while (norm(v_k1 - v_k) / norm(v_k)) > tol
    v_k = v_k1;
    w = A*v_k1;
    v_k1 = w/max(w);
end
eigvec = v_k1;
eigval = max(w);
end