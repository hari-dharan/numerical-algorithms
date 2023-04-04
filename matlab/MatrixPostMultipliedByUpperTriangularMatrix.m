function C = MatrixPostMultipliedByUpperTriangularMatrix(A, B)
%MATRIXPOSTMULTIPLYUPPERTRIANGULARMATRIX Post multiply with an upper triangular matrix.
%   C = A * B where C is (mxn), A is (mxn) and B is (nxn) and B is upper triangular.
%   Redundant multiplications with zero elements in B are avoided.
%   Number of multiplications: mn(n+1)/2
%   Number of additions: mn(n-1)/2
m = size(A,1);
n = size(A,2);
C = zeros(m,n);
for i=1:m
    for j=1:n
        C(i,j) = A(i,1)*B(1,j);
        for k=2:j
            C(i,j) = C(i,j) + A(i,k)*B(k,j);
        end
    end
end
end