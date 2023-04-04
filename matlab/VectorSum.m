function sum = VectorSum(x)
%VECTORSUM sums the elements of a vector x.
%   Uses naive for loop to sum elements iteratively.
sum = 0;
n = length(x);
for i=1:n
    sum = sum + x(i);
end
end