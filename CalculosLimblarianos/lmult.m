function S = lmult(A)
% Super-operator representing A being left-multiplied:
% S = ()*A

% Given a (N x M) matrix A this function constructs a new (N^2 x NM) matrix
% S. Given another (K x N) matrix B if this is flatterned into a (KN x 1)
% column vector b then S*b returns the flatterned (KM x 1) column vector a
% that would result from left-multipling B by A as B*A.

% Get dimension of operator:
d = size(A,2);

S = kron(sparse(A.'), speye(d));