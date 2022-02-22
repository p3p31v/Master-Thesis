%% diss.m by Mark Mitchison
% Computes a matrix representation of a dissipator

function A = diss(L,B)

A = L*B*L' - 1/2*L'*L*B - 1/2*B*L'*L;

end