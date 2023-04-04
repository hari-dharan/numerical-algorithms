function sum = VectorSumKahan(x)
%VECTORSUMKAHAN sums the elements of a vector x using Kahan's Algorithm.
%   Naive sum loses accuracy with a large vector x. Kahan's Algorithm 
%   uses a correction variable to keep track of the error during each 
%   summation in the loop. The effect is increased accuracy but slower 
%   speed due to more operations.
n = length(x);
sum = 0;
c = 0;
for i=1:n
    y = x(i) - c;
    t = sum + y;
    c = (t-sum) - y;
    sum = t;
end
end