%% lindDeriv.m 
% Calculates the action of a Lindblad generator on a quantum state

function rhoOut = lindDeriv(H, L, rhoIn)

% Number of Lindblad operators
nLind = length(L);

% Coherent part
rhoOut = -1i*commutator(H,rhoIn);

% Incoherent part
for j = 1:nLind
    rhoOut = rhoOut + diss(L{j},rhoIn);
end

end