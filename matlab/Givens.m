function y = Givens(x, i, j)
%GIVENS Use Givens rotation to eliminate one entry of a vector.
%   x: input column vector.
%   i: coordinate of the entry used to eliminate
%   j: coordinate of the entry to eliminate
%   The Givens rotation matrix is not explicitly formed. 
%   The output is the vector with the j-th entry eliminated.
a = x(i);
b = x(j);

if b == 0
    c = 1;
    s = 0;
else
    if abs(b) > abs(a)
        r = a / b;
        s = 1 / sqrt(1 + r^2);
        c = s * r;
    else
        r = b / a;
        c = 1 / sqrt(1 + r^2);
        s = c * r;
    end
end

y = x;
y(i) = c*a + s*b;
y(j) = -s*a + c*b;

end