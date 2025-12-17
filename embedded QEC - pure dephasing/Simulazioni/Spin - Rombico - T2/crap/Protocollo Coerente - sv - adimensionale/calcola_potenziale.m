function V = calcola_potenziale(sys, imp, t)
%   Funzione per il calcolo della forma opportuna del potenziale in
%   interaction picturedovuto all'impulso passato

%%% Per eseguire le simulazioni in interaction picture
%   Calcoliamo l'operatore da sch ad int
    U = diag(exp(1i * mod(diag(sys.H0) * t, 2*pi)));

% % % %%% Per eseguire le simulazini in schrodinger picture
% % %     U = eye(sys.dim);

%   Parte risonante lungo z
    V = (imp.Bz * sys.gz * sys.mu_b * cos(mod(imp.w * t + imp.fase, 2*pi)));
end