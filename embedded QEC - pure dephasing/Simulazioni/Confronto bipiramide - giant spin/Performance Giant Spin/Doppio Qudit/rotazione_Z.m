function R = rotazione_Z(phi)
%   Funzione che ritorna la rotazione piana di parametri theta e beta
    R = [1 0; 0 exp(1i * phi)];
end