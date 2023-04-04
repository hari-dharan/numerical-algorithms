function x = ForwardSubstitution(A)
%FORWARDSUBSTITUTION Forward substitution of an augmented lower triangular matrix.
n = size(A,1);
x = zeros(n,1);
x(1) = A(1,n+1)/A(1,1);
for i=2:n
    x(i) = A(i,n+1);
    for j=1:i-1
        x(i) = x(i) - A(i,j)*x(j);
    end
    x(i) = x(i)/A(i,i);
end
end