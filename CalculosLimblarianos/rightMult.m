% rightMult.m by Mark T. Mitchison (C) Imperial College London 2014

% Takes a matrix A and returns a superoperator representing the action of
% multiplying from the right on a flattened density operator. The density operator
% should be represented by the column vector reshape(rho,dim^2,1), where dim is the 
% dimension of the Hilbert space.

function super = rightMult(A)

% Size of the matrix A
dimension = size(A,1);

% Generate superoperator
super = kron(A.',speye(dimension));
