function x = BackSubstitution(A)
%BACKSUBSTITUTION Back substitution of an augmented upper triangular matrix.
n = size(A,1);
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