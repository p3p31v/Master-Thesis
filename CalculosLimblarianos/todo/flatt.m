% flatt.m by Mark T. Mitchison (c) Imperial College London 2014
% Reshapes a square matrix into a column vector.

function vector = flatt(matrix) 

vecLength = prod(size(matrix)) % Length of the new vector

vector = reshape(matrix,vecLength,1) % Reshape into vector

endfunction