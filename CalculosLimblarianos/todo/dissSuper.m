%% dissSuper.m by Mark Mitchison (C) Imperial College London 2015
% Generates a sparse matrix representation of a dissipator superoperator with a 
% Lindblad operator L

function diss = dissSuper(L)

    diss = leftMult(L)*rightMult(L') - 1/2*leftMult(L'*L) - 1/2*rightMult(L'*L);
    
end