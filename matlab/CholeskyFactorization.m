function L = CholeskyFactorization(A)
%CHOLESKYFACTORIZATION Cholesky Factorization of a symmetric positive definite matrix A.
%   Cholesky Factorization yields A = LL' where L' is the transpose of L.
%   This function returns L.
%   Number of square roots: n
%   Number of divisions: n(n-1)/2
%   Number of multiplications: n^2(n-1) - (2n(n)(n-1))/2 + (n(n-1)(2n-1))/6
%   Number of subtractions: Same as number of multiplications
n = size(A,1);
for i=1:n-1
    A(i,i) = sqrt(A(i,i));
    for j=i+1:n
        A(j,i) = A(j,i)/A(i,i);
    end
    for j=i+1:n
        for k=i+1:n
            A(k,j) = A(k,j) - A(k,i)*A(j,i);
        end
    end
end
A(n,n) = sqrt(A(n,n));
L = tril(A);
end