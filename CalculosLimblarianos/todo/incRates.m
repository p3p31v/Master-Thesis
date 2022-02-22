
function gamma = incRates(E,wMax,alpha,s,kDiag,beta,mu,fermions);

% Generate incoherent rates
gamma = zeros(3,3);
for j = 1:3 % Loop over sites
    % Calculate rates for fermionic bath
    if(fermions(j))
        gamma(j,1) = spectPow(E(j),wMax(j),alpha(j),s(j))*fermiDist(E(j),beta(j),mu(j)); % Gain rate
        gamma(j,2) = spectPow(E(j),wMax(j),alpha(j),s(j))*(1-fermiDist(E(j),beta(j),mu(j))); % Loss rate
        gamma(j,3) = 0; % Dephasing rate
    % Calculate rates for bosonic bath
    else 
        gamma(j,1) = spectPow(E(j),wMax(j),alpha(j),s(j))*boseDist(E(j),beta(j)); % Gain rate
        gamma(j,2) = spectPow(E(j),wMax(j),alpha(j),s(j))*(1+boseDist(E(j),beta(j))); % Loss rate
        if s == 1
            gamma(j,3) = 2*kDiag(j)*alpha(j)/beta(j); % Dephasing rate only finite for an Ohmic bosonic bath
        else 
            gamma(j,3) = 0;
        end
    end
end
end

% Power law spectral density with exponential cutoff
function y = spectPow(w,wc,alpha,s)
y = alpha*w^s*exp(-w/wc);
end

% Bose-Einstein distribution
function y = boseDist(w,beta)
y = (exp(beta*w) - 1)^-1;
end

% Fermi-Dirac distribution
function y = fermiDist(w,beta,mu)
y = (exp(beta*(w-mu))+1)^-1;
end

