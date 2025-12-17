function dy = forzante_lindblad(t,y)
%   Funzione per l'implemntazione del forzante dell'equazione di Lindblad
%   per l'evoluzione del sistema in schrodinger picture

%   Variabile globale che contiene tutti i dati relativi alla simulazione
%   da eseguire
    global param;
    
%   Dimensione del sistema
    dim = param.sys.dim;

%   Recuperiamo la corretta forma per la rho
    rho = reshape(y, dim, dim);

%   Recuperiamo il numero di impulsi da simulare in parallelo
    n_imp = length(param.imp);

%   Inizializzazione Variabili
    H_tot = zeros(dim,dim);
    drho  = zeros(dim,dim);

% % % %%% Per simuazioni in schrodinger
% % % %   Includiamo l'Hamiltoniana
% % %     H_tot = param.sys.H0;

%   Calcoliamo la parte coerente dell'evoluzione
    for i = 1:n_imp
%       Calcoliamo la forma opportuna del potenziale dovuto all'impulso
%       i-esimo in interaction picture
        V_I = calcola_potenziale(param.sys, param.imp{i}, t);

%       Aggiungiamo il contributo del potenziale all'evoluzione del sistema
        H_tot = H_tot + V_I;
    end

%   Aggiungiamo il contribu

%   Aggiungiamo il contributo del termine statico
    H_tot = H_tot + calcola_potenziale_statico(param, t);

%   Includiamo la parte coerente nell'evoluzione del sistema
    drho = drho + 1i * (rho * H_tot - H_tot * rho);

%   Includiamo la parte incoerente nell'evoluzione del sistema
%     drho = drho + pure_dephasing(param.sys, rho);
    
%   Ripristiniamo la corretta forma per l'output della funzione
    dy = reshape(drho, dim^2, 1);
end