function R = rotazione_piana(theta, beta)
%   Funzione che ritorna la rotazione piana di parametri theta e beta
    R = [cos(theta) -sin(theta)*exp(-1i*beta); sin(theta)*exp(1i*beta) cos(theta)];
end