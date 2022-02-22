%% bosonOpsInhom.m by Mark Mitchison (c) Imperial College London 2015
% Generates a set of sparse bosonic operators on L sites with varying 
% occupation cutoffs.

function [B, Bd, Num, Id] = bosonOpsInhom(L, nMax)

% If system is homogeneous make nMax into a vector
if length(nMax) == 1
    nMax = nMax*ones(1,L);
elseif length(nMax) ~= L
    error('nMax must have length L!')
end

% Create boson operators for each site
for l = 1:L
    
    % Dimension of the local Hilbert space
    dim = nMax(l)+1; 

    % Vector of off-diagonal elements of the annihilation operator
    sqrtNum = sqrt(0:nMax(l)+1)';    
    
    % Define single-site operators
    b = spdiags(sqrtNum, 1, dim, dim);
    n = spdiags((0:nMax(l))', 0, dim, dim);
    
    % Create the identity operator corresponding to sites on the left
    Il = 1;
    for m = 1:l-1
        Il = kron(Il,speye(nMax(m)+1));
    end
    % Create the identity operator corresponding to sites on the right
    Ir = 1;
    for m = l+1:L
        Ir = kron(Ir,speye(nMax(m)+1));
    end
    
    % Tensor single-site operators into full system
    B{l} = kron(kron(Il,b),Ir);
    Num{l} = kron(kron(Il,n),Ir);
    Bd{l} = B{l}';
end

% Create full identity
Id = kron(Il, speye(nMax(L)+1));
