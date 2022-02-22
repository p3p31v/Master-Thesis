%% adjDiss.m by Mark Mitchison (C) Imperial College London 2015
% Generates a sparse operator representation of an adjoint dissipator with 
% Lindblad operator L

function diss = adjDiss(L)

    diss = leftMult(L')*rightMult(L) - 1/2*leftMult(L'*L) - 1/2*rightMult(L'*L);
    
end