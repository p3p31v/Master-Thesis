%% bosonOps.m by Mark Mitchison (C) Imperial College 2014
% Creates a set of L sparse bosonic creation and annihilation operators,
% admitting a maximum of nMax bosons per mode. Also creates the number 
% operator for convenience.

function [B, Bd, Num, Id] = bosonOps(L,nMax)

% Dimension of the local Hilbert space
dim = nMax+1; 

% Vector of off-diagonal elements of the annihilation operator
sqrtNum = sqrt(0:nMax+1)';

% Create one-site boson operators
b = spdiags(sqrtNum, 1, dim, dim);
num = spdiags((0:nMax)', 0, dim, dim);

% Tensor together to give the full set
for j = 1:L
    B{j} = kron(speye(dim^(j-1)), kron(b, speye(dim^(L-j))));
    Bd{j} = B{j}';
    Num{j} = kron(speye(dim^(j-1)), kron(num, speye(dim^(L-j))));
end

Id = speye(dim^L);