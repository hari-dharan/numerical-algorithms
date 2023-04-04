function x = SolveByGaussianElimination(A, b)
%SOLVEBYGAUSSIANELIMINATION Solve linear systems of equations Ax=b.
%   Solved using Gaussian Elimination then Backward Substitution.
%   A is (nxn) and b is (nx1).
%   The algorithm leaves eliminated entries as the original entries.
%   Number of divisions: n(n+1)/2
%   Number of multiplications: n(n-1)(2n+5)/6
%   Number of subtractions: n(n-1)(2n+5)/6

% Create augmented matrix
A = [A b];

% Gaussian Elimination
n = size(A,1);
for i=1:n-1
    for j=i+1:n
        m = A(j,i)/A(i,i);
        for k=i:n+1
            A(j,k) = A(j,k) - m*A(i,k);
        end
    end
end

% Backward Substitution
x = zeros(n,1);
x(n) = A(n,n+1)/A(n,n);
for i=n-1:-1:1
    x(i) = A(i,n+1);
    for j=i+1:n
        x(i) = x(i) - A(i,j)*x(j);
    end
    x(i) = x(i)/A(i,i);
end

end