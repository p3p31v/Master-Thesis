%% lindbladGen.m 
% Calculates the action of a Lindblad generator on a quantum state

function rhoOut = lindbladGen(H, L, rhoIn)

% Number of Lindblad operators
nLind = length(L);

% Coherent part
rhoOut = -1i*commutator(H,rhoIn);

% Incoherent part
for j = 1:nLind
    rhoOut = rhoOut + diss(L,rhoIn);
end

end