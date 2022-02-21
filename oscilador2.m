% Simulation parameters
ss = 0;                     % Set flag to 1 to find the stationary state
tTot = 4000;                 % Total simulation time
dt = 0.1;                           % Numerical time step
nSteps = ceil(tTot/dt);           % Number of time steps
tol = 1e-10 ;                    % Numerical tolerance
time = linspace(0,tTot,nSteps+1);  % Time axis

%% Generators and states




psi0=kron([1,0,0],kron([1,0],[0,1]));

%psi0=kron([1,0,0],kron(kron([1,0],[0,1])+kron([0,1],[1,0]),[1,0,0]));

%psi0=(1/sqrt(2))*(kron([1,0,0],kron(kron([1,0],[0,1])+kron([0,1],[1,0]),[1,0,0])));

%psi0=(kron([1,0,0],0.9239*kron([1,0],kron([0,1],[1,0,0])))+0.3827*kron([1,0,0],kron([0,1],kron([1,0],[1,0,0]))));

%psi0=(1/sqrt(2))*kron([1,0,0],kron([1,0],[0,1]))+(1/sqrt(2))*kron([1,0,0],kron([0,1],[1,0]));




rho0=psi0'*psi0;

size(rho0);

rho0=rho0/trace(rho0);


[B, Bd, Num, Id] = bosonOpsInhom(3,[2,1,1]);

Sx1 = Bd{1}+B{1};
Sy1 = 1i*(B{1}-Bd{1});
Sz1 = 2*Num{1}-Id;
Sx2 = B{2}+Bd{2};
Sy2= 1i*(B{2}-Bd{2});
Sz2=2*Num{2}-Id;
Sx3=Bd{3}+B{3};
Sy3=1i*(B{3}-Bd{3});
Sz3=2*Num{3}-Id;


j12=sqrt(2);
j23=100;

j13=0;

w1=200;


delta2=0;
delta3=0;


H = delta2*Sz2/2+delta3*Sz3/2+j12*(B{1}+Bd{1})*(Sz2+Id)+j23*(Bd{2}*B{3}+B{2}*Bd{3})+w1*Num{1}+j13*(B{1}+Bd{1})*(Sz3+Id);

%nosc
[Bo, Bdo, Numo, Ido] = bosonOpsInhom(1,2);

Sx1o = Bdo{1}+Bo{1};
Sy1o = 1i*(Bo{1}-Bdo{1});
Sz1o = 2*Numo{1}-Ido;


delta1osc=delta2;
delta2osc=delta3;
j12osc=j23;


Hosc =w1*Numo{1};


Hosc=full(Hosc);
lambdao=eig(Hosc);
[u,lambdao]=eig(Hosc);
I=eye(3);
In=eye(2)
nosc1=kron(u(:,1)*u(:,1)',kron(In,In));
nosc2=kron(u(:,2)*u(:,2)',kron(In,In));
nosc3=kron(u(:,3)*u(:,3)',kron(In,In));
coheosc=kron(u(:,1)*u(:,2)',kron(In,In));




%hasta aquÃ­



%nexc

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
I=eye(3);
nexc1=kron(I,v(:,2)*v(:,2)');
nexc2=kron(I,v(:,3)*v(:,3)');
cohe12=kron(I,v(:,1)*v(:,4)');
cohe1g=kron(I,v(:,2)*v(:,4)');
cohe2g=kron(I,v(:,3)*v(:,4)');


%hasta aquÃ­



% Coherent Liouvillian
L = -1i*lmult(H) + 1i*rmult(H);


% Lindblad operators
nLind = 0 ;

lambda=0 ;                                % Number of Lindblad operators
if lambda ~= 0

	lind{1}=0.03*Sz2;
	lind{2}=0.03*Sz3;   
end

% Add dissipators
for j = 1:nLind
    L = L + dissSuper(lind{j});
end

%% Stationary state






    
rho = rho0;

for ind = 1:nSteps
    % Compute observables

    population(ind)=trace(rho*Num{3});
	
populationosc1(ind)=trace(nosc1*rho);

populationosc2(ind)=trace(nosc2*rho);

populationosc3(ind)=trace(nosc3*rho);



    populationexc2(ind)=trace(nexc2*rho);

    populationexc1(ind)=trace(nexc1*rho);

    coherence12(ind)=real(trace(cohe12*rho));

    coherence1g(ind)=real(trace(cohe1g*rho));

    coherence2g(ind)=real(trace(cohe2g*rho));
    
    coherenceosc(ind)=real(trace(coheosc*rho));

    % Evolve forward one time step

    rho = evolveOp(dt, L, rho, tol);
    

    
end





save variable5.out populationosc1 -ASCII

save variable6.out populationosc2 -ASCII

save variable7.out populationosc3 -ASCII

save variable8.out coherenceosc -ASCII

save variable.out populationexc1 -ASCII

save variable2.out populationexc2 -ASCII

save variable3.out coherence12 -ASCII 

save variable4.out coherence1g -ASCII 

