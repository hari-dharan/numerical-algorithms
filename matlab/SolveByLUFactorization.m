function x = SolveByLUFactorization(L, U, b)
%SOLVEBYLUFACTORIZATION Summary of this function goes here
%   Detailed explanation goes here
n = size(L,1);

% Solve Lb' = b using forward substitution
L = [L b]; % augmented matrix L
b_ = zeros(n,1); % b' is variable b_
b_(1) = L(1,n+1)/L(1,1);
for i=2:n
    b_(i) = L(i,n+1);
    for j=1:i-1
        b_(i) = b_(i) - L(i,j)*b_(j);
    end
    b_(i) = b_(i)/L(i,i);
end

% Solve Ux = b' using backward substitution
U = [U b_]; % augmented matrix U
x = zeros(n,1);
x(n) = U(n,n+1)/U(n,n);
for i=n-1:-1:1
    x(i) = U(i,n+1);
    for j=i+1:n
        x(i) = x(i) - U(i,j)*x(j);
    end
    x(i) = x(i)/U(i,i);
end

end