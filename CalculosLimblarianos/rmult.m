function S = rmult(A)
% Super-operator representing A being right-multiplied:
% S  = A*()

% Given a (N x M) matrix A this function constructs a new (N^2 x NM) matrix
% S. Given another (M x K) matrix B if this is flatterned into a (MK x 1)
% column vector b then S*b returns the flatterned (NK x 1) column vector a
% that would result from right-multipling B by A as A*B.

% Get dimension of operator:
d = size(A,2);

S = kron(speye(d), sparse(A));