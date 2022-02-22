%% dissipator.m by Mark Mitchison (C) Imperial College London 2015
% Generates a sparse operator representation of a dissipator with a rate
% Gamma and a Lindblad operator L

function diss = dissipator(L)

    diss = leftMult(L)*rightMult(L') - 1/2*leftMult(L'*L) - 1/2*rightMult(L'*L);
    
end