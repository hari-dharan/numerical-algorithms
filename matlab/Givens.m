function G = Givens(x, i, j)
%GIVENS Return Givens rotation matrix to eliminate one entry of a vector.
%   x: input column vector.
%   i: coordinate of the entry used to eliminate
%   j: coordinate of the entry to eliminate 
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

G = eye(length(x));
G(i,i) = c;
G(j,j) = c;
G(i,j) = s;
G(j,i) = -s;

end