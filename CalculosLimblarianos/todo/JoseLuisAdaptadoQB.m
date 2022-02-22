






% Simulation parameters
ss = 0;                             % Set flag to 1 to find the stationary state
tTot = 3000;                        % Total simulation time
dt = 1;                             % Numerical time step
nSteps = ceil(tTot/dt);             % Number of time steps
tol = 1e-10;                         % Numerical tolerance
time = linspace(0,tTot,nSteps+1);   % Time axis

%% Generators and states

% Qubit thermal state
rhoQ = [1 0; 0 exp(-beta)];
rhoQ = rhoQ/trace(rhoQ);

% Coherent initial battery state
rhoB = zeros(1,d)';
for j = 1:d
    rhoB(j) = aBatt^(j-1);
end

rhoB = rhoB*rhoB';
rhoB = rhoB/trace(rhoB);

% Initial state
rho0 = kron(rhoQ,rhoB);


[B, Bd, Num, Id] = bosonOps(2,[1,d-1]);

Sx = Bd{1}+B{1};
Sy = 1i*(B{1}-Bd{1});
Sz = 2*Num{1}-Id;
X = B{2}+Bd{2};
Y = 1i*(B{2}-Bd{2});


for jj=1:200

H = Omega/2*B{1} + g*B{1}*Bd{2} - cos(Ang)*f*Bd{1}*B{1} ; 

H = H + H';

% Coherent Liouvillian
L = -1i*leftMult(H) + 1i*rightMult(H);

% Lindblad operators
nLind = 1;                                      % Number of Lindblad operators
if lambda ~= 0

    lind{1} = Sz; % Measurement and
   
end

% Add dissipators
for j = 1:nLind
    L = L + dissSuper(lind{j});
end

%% Stationary state

if ss   % Compute stationary state

    % Find zero eigenvector
    [rhoSS,eVal] = eigs(L,1,'sm');
    rhoSS = unflatt(rhoSS);
    rhoSS=rhoSS/trace(rhoSS);
    
    % Compute observables
    sigxSS = trace(rhoSS*Sx);
    sigxSS = trace(rhoSS*Sy);
    sigxSS = trace(rhoSS*Sz);
    NbattSS = trace(rhoSS*Num{2});
    
else 
    
rho = rho0;

% Initial observables
sigx = zeros(1,nSteps+1); sigy = zeros(1,nSteps+1); sigz = zeros(1,nSteps+1); Nbatt = zeros(1,nSteps+1); 

sigx(1) = trace(rho*Sx);
sigy(1) = trace(rho*Sy);
sigz(1) = trace(rho*Sz);

for ind = 1:nSteps
    
    % Evolve forward one time step
    rho = evolveOp(dt, L, rho, tol);
   
    % Compute observables
    sigx(ind+1) = trace(rho*Sx);
    sigy(ind+1) = trace(rho*Sy);
    sigz(ind+1) = trace(rho*Sz); 
    
end


end % Este es el bucle que esta escribiendo las distintas f


%% Plots

% Battery population
%figure
hold on
%save (time, Nbatt, time, ergB);
title('med sx lamb = 1 , eta = 1, dt = 1, Bat = 30, T = 1, Gam = 1 TW-B copling')
save ('Trisg10Lam1Eta1G1Bat30T1', 'time', 'Nbatt', 'ergB');
%plot(time, Nbatt, time, ex, time, wy, time, ergB)

legend('g=2')

plot(time, Nbatt)


legend('g=2')

% Qubit observables
%%figure % Se necesita para que abra otra figura
hold on


end