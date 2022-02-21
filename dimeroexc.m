% Simulation parameters
ss = 0 ;                         % Set flag to 1 to find the stationary state
tTot = 1 ;                     % Total simulation time
dt = 0.001    ;                       % Numerical time step
nSteps = ceil(tTot/dt)  ;          % Number of time steps
tol = 1e-10          ;            % Numerical tolerance
time = linspace(0,tTot,nSteps+1) ; % Time axis



% Initial state

psi0=kron([0,1],[1,0]);

rho0=psi0'*psi0;

size(rho0);

rho0=rho0/trace(rho0);



[B, Bd, Num, Id] = bosonOpsInhom(2,[1,1]);

Sx1 = Bd{1}+B{1};
Sy1 = 1i*(B{1}-Bd{1});
Sz1 = 2*Num{1}-Id;
Sx2 = B{2}+Bd{2};
Sy2= 1i*(B{2}-Bd{2});
Sz2=2*Num{2}-Id;


delta1=20;
delta2=20;
j12=8;


H = delta1*Sz1/2+delta2*Sz2/2+j12*(Bd{1}*B{2}+B{1}*Bd{2});


[Be, Bde, Nume, Ide] = bosonOpsInhom(2,[1,1]);

Sx1e = Bde{1}+Be{1};
Sy1e = 1i*(Be{1}-Bde{1});
Sz1e = 2*Nume{1}-Ide;
Sx2e = Be{2}+Bde{2};
Sy2e= 1i*(Be{2}-Bde{2});
Sz2e=2*Nume{2}-Ide;


delta1exc=delta2;
delta2exc=delta3;
j12exc=j23;


Hexc = delta1exc*Sz1e/2+delta2exc*Sz2e/2+j12exc*(Bde{1}*Be{2}+Be{1}*Bde{2});



Hexc=full(Hexc);
lambdae=eig(Hexc);
[v,lambdae]=eig(Hexc);
nexc1=kron(v(:,2)*v(:,2)');
nexc2=kron(v(:,3)*v(:,3)');
cohe12=kron(v(:,2)*v(:,3)');



% Coherent Liouvillian
L = -1i*lmult(H) + 1i*rmult(H);



% Lindblad operators
nLind = 0      ;                         % Number of Lindblad operators

lambda=0;
if lambda ~= 0

    lind{2} = 3*Bd{1};
	lind{1}=0.9*Sz1;
	lind{3}= 0.3*Sz2;
	lind{4}=0.4*Bd{2};
   
end

% Add dissipators
for j = 1:nLind
    L = L + dissSuper(lind{j});
end



rho=rho0;
for ind = 1:nSteps
    
    % Compute observables
    nuuu(ind) = trace(rho*Num{1});
nuuuu(ind) = trace(rho*Num{2});


    populationexc2(ind)=trace(nexc2*rho);

    populationexc1(ind)=trace(nexc1*rho);

    coherence12(ind)=real(trace(cohe12*rho));


    % Evolve forward one time step
    rho = evolveOp(dt, L, rho, tol);

    
    
end

save variable.out nuuu -ASCII
save variable2.out nuuuu -ASCII

save variable.out populationexc1 -ASCII

save variable2.out populationexc2 -ASCII

save variable3.out coherence12 -ASCII 

