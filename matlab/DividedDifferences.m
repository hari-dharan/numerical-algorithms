function b = DividedDifferences(x, y)
%DIVIDEDDIFFERENCES Compute the coefficients of the Newton Basis.
n = length(x);
b = y;
for j = 2:n
    for k = n:-1:(j+1)
        b(k) = (b(k) - b(k-1)) / (x(k) - x(k-j));
    end
end
end