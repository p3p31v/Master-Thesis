%% evolveOp.m by Mark Mitchison 
% Evolves an operator forward in time by exponentiatin a generator using 
% the Lanczos technique

function [opOut] = evolveOp(dt, G, opIn, tol)

opOut = flatt(opIn);

opOut = expv(dt, G, opOut, tol);

opOut = unflatt(opOut);

end
