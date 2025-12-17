function dy = forzante_lindblad(t,y)
%   Funzione per l'implemntazione del forzante dell'equazione di Lindblad
%   per l'evoluzione del sistema in schrodinger picture

%   Variabile globale che contiene tutti i dati relativi alla simulazione
%   da eseguire
    global param;
    
%   Dimensione del sistema
    dim = param.sys.dim;

%   Recuperiamo il numero di impulsi da simulare in parallelo
    n_imp = length(param.imp);

%   Inizializzazione Variabili
    H_tot = zeros(dim,dim);
    V_I   = 0;

% % % %%% Per simuazioni in schrodinger
% % % %   Includiamo l'Hamiltoniana
% % %     H_tot = param.sys.H0;

%   Calcoliamo la parte coerente dell'evoluzione
    for i = 1:n_imp
%       Calcoliamo la forma opportuna del potenziale dovuto all'impulso
%       i-esimo in interaction picture
        V_I = V_I + calcola_potenziale(param.sys, param.imp{i}, t);
    end

%   Rimozione parziale leackage
    U = diag(exp(1i * mod(diag(param.sys.H0) * t, 2*pi)));
    H_tot = V_I * (U * diag(diag(param.sys.Mz)) * U');

%   Aggiungiamo il contributo del termine statico
    H_tot = H_tot + calcola_potenziale_statico(param, t);

%   Includiamo la parte coerente nell'evoluzione del sistema
    dy = - 1i * (H_tot) * y;
end