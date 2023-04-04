function C = MatrixMultiplication(A, B)
%MATRIXMULTIPLICATION Product of two matrices under the usual matrix multiplication.
%   C = A * B where C is (mxp), A is (mxn) and B is (nxp).
%   Naive algorithm of 3 for loops to compute each entry in output matrix.
m = size(A, 1);
n = size(A, 2);
p = size(B, 2);
C = zeros(m, p);
for i=1:m
    for j=1:p
        C(i,j) = A(i,1)*B(1,j);
        for k=2:n
            C(i,j) = C(i,j) + A(i,k)*B(k,j);
        end
    end
end
end