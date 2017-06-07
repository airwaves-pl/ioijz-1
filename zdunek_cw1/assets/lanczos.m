function T = lanczos(A,v,k)
% LANCZOS A simple implementation of the Lanczos process using matrixop
%
% T = lanczos(A,v,k) returns the tridiagonal matrix T of coefficients for
% the Lanczos process after k steps, starting with the vector v.  The
% matrix A must be a matrixop class or a Matlab matrix.  The operator
% matrix A must be symmetric.
%
% See Golub and van Loan, Matrix Computations, 3rd Edition for more
% information.

% initialization
v = v./norm(v);
r = v;
b1 = 1;
v1 = 0;
T = zeros(k,k);

% computation
i = 0;
while (i < k-1)
    v = r/b1;
    i = i+1;
    p = A*v;
    T(i,i) = v'*p;
    r = p - T(i,i)*v - b1*v1;
    b1 = norm(r);
    T(i,i+1) = b1;
    v1 = v;
end
% finish the final step
v = r/b1;
i = i+1;
p = A*v;
T(i,i) = v'*p;

% set the lower diagonal
T = T + diag(diag(T,1),-1);