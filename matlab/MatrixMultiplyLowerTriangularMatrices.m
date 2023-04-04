function C = MatrixMultiplyLowerTriangularMatrices(A, B)
%MATRIXMULTIPLYLOWERTRIANGULARMATRICES Product of two lower triangular matrices 
% under the usual matrix multiplication.
%   C = A * B where C is (nxn), A is (nxn) and B is (nxn) and A, B are lower triangular.
%   Redundant multiplications with zero elements in A and B are avoided.
%   Number of multiplications: n(n+1)(n+2)/6
%   Number of additions: (n-1)n(n+1)/6
n = size(A,1);
C = zeros(n,n);
for i=1:n
    for j=1:i
        C(i,j) = A(i,j)*B(j,j);
        for k=j+1:i
            C(i,j) = C(i,j) + A(i,k)*B(k,j);
        end
    end
end
end