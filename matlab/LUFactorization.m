function [L,U] = LUFactorization(A)
%LUFACTORIZATION Return LU factors of given (nxn) square matrix A.
%   Same algorithm as Gaussian Elimination except the multipliers are
%   stored.
n = size(A,1);
for i=1:n-1
    for j=i+1:n
        A(j,i) = A(j,i)/A(i,i);
        for k=i+1:n
            A(j,k) = A(j,k) - A(j,i)*A(i,k);
        end
    end
end
L = tril(A,-1);
U = triu(A);
end