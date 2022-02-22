function n  = BEdist(E,T)
% BEdist calculates the Bose-Einstein distribtion for energy E and
% temperature T

n = 1/(exp(E/T) - 1);

end

