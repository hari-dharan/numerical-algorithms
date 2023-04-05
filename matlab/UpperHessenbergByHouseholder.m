function hess = UpperHessenbergByHouseholder(A)
%UPPERHESSENBERGBYHOUSEHOLDER Create an upper hessenberg matrix.
n = size(A,1);
hess = A(2:n, 1:n-1);
n_hess = size(hess,1);
HH1 = Householder(hess(:,1), 0);
hess = HH1*hess;
for i=2:n_hess
    HH = Householder(hess(i:n_hess, i), i-1);
    hess = HH*hess;
end
hess = [A(1,:); hess, A(2:n,n)];
end