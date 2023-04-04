function A = GaussianElimination(A)
%GAUSSIANELIMINATION Gaussian Elimination of an augmented matrix A, (nx(n+1)).
%   Usual Gaussian Elimination of an augmented matrix.
n = size(A,1);
for i=1:n-1
    for j=i+1:n
        m = A(j,i)/A(i,i);
        for k=i:n+1
            A(j,k) = A(j,k) - m*A(i,k);
        end
    end
end
end