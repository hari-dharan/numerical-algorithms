function [eigvec, eigval, n_iter] = RayleighQuotientIteration(A, v, tol)
% v_k is v at k-th iteration and v_k1 is v at k+1-th iteration
v_k = zeros(size(v));
v_k1 = v;
n_iter = 0;
while (norm(v_k1 - v_k) / norm(v_k)) > tol
    v_k = v_k1;
    sigma = (v_k1'*A*v_k1) / (v_k1'*v_k1);
    w = linsolve(A-sigma*eye(size(A)), v_k1);
    v_k1 = w/max(w);
    n_iter = n_iter + 1;
end
sigma = (v_k1'*A*v_k1) / (v_k1'*v_k1);
eigvec = v_k1;
eigval = sigma;