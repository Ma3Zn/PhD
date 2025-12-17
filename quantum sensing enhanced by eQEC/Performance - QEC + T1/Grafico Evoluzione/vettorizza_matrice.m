function [v] = vettorizza_matrice(rho)
    n = length(rho);
    v = reshape(rho.', n^2, 1);
end