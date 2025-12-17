function [rho] = devettorizza_matrice(v)
    n = length(v);
    n = n^0.5;

    rho = reshape(v, n, n).';
end