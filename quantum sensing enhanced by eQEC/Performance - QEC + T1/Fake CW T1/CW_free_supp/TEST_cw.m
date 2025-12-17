clear
clc
tic

%   script per la ricerca di codewords numeriche per la protezione di
%   sistemi con interazioni in competizione solo per T2 (Sz). Usiamo un
%   encoding reale

dim = 4;

E = calcola_operatori_errore(dim, 5e-6);

%   Generazione parametri errore
global ERR;

%   Calcoliamo le dimensioni del sistema
ERR.dim  = dim;
ERR.n    = ERR.dim / 2;

ERR.E = cell(1,ERR.n);

%   Inseriamo gli operatori di Krauss che vogliamo correggere
for j = 1:ERR.n
    ERR.E{j} = E{j};
end

%   Inseriamo lo spacing per le codewords. sempre 1 per T2
ERR.spacing = 1;

% % % ERR.supp0 = 1:2:ERR.dim;
% % % ERR.supp1 = ERR.dim + [2:2:ERR.dim];

%   Creazione handle per la funzione da ottimizzare e pe ivincoli non
%   lineari
nonlcon = @nonlinear_constraint;

%   Generiamo upeer e lower bound per i valori delle variabili di encoding
lb = -1 * ones(2*ERR.dim,1);
ub =      ones(2*ERR.dim,1);

%   ATTENZIONE: Possibile sia uno dei punti più delicati di tutto
%   l'algoritmo, i.e., il punto iniziale
x_init = genera_condizione_iniziale_STD(ERR.dim);
% x_init = genera_condizione_iniziale_MOD(ERR.dim);

%   Settiamo le opzioni di ricerca per il solver
options = struct();
options.maxfun = 5e2*length(x_init);
% options.rhoend = eps;
% options.rhobeg = 1e-3;

%   Eseguiamo l'otimizzazione utilizzando l'interfaccia MATLAB di PRIMA,
%   pacchetto moderno di ottimizzatori (Powell's algorithms) scritto in
%   fortran
[x_opt, ~, exitflag, output] = prima([], x_init, [], [], [], [], lb, ub, nonlcon, options);

ew0 = x_opt(1:ERR.dim);
ew1 = x_opt(ERR.dim+1:end);

[a,b] = nonlcon(x_opt);
[c,d] = nonlcon(x_init);

exitflag;

output;

toc