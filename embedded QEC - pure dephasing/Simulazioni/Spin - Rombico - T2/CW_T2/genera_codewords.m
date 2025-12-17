function [ew0, ew1, pass] = genera_codewords(dim, tol)
%   Funzione che data in input la dimensione del sistema ed il tempo di
%   ottimizzazione genera le codewords per il sistema tramite un'iterazione
%   del processo di ottimizzazione. Iterazione che si spera porti ad un via
%   via miglioramento delle codewords.

%   Controlliamo la consistenza della tolleranza
    if tol < 1e-15
        tol = 1e-15;
    end
    
%   Generiamo le variabili di output
    pass = false;
    cw   = zeros(dim,2);

%   Settiamo gli operatori d'erore
    E = operatori_errori(dim);

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

%   Supporti delle codifiche
    global supp0;
    global supp1;

%   Generiamo la condizione iniziale
    x_init    = genera_condizione_iniziale(dim);

%   Creazione handle per la funzione da ottimizzare e per i vincoli non
%   lineari
    nonlcon = @nonlinear_constraint;

%   Generiamo upper e lower bound per i valori delle variabili di encoding
    lb = -1 * ones(ERR.dim,1);
    ub =      ones(ERR.dim,1);

%   Settiamo le opzioni di ricerca per il solver
    options = struct();
    options.maxfun = 5e2*length(x_init);

%   Eseguiamo il primo step di ottimizzazione
    [x_opt, ~, ~, output] = prima([], x_init, [], [], [], [], lb, ub, nonlcon, options);

%   Normalizziamo x_opt
    x_opt = pulisci_codewords(x_opt);

%   Contatore per il numero di iterazioni del ciclo di ottimizzazione
    it_opt = 1;

%   Iteriamo l'ottimizzazine
    while (output.constrviolation > tol && it_opt < 10)

%       Aggiorniamo la condizione iniziale con le codewords trovate al
%       passo precedente
        x_init = x_opt;

%       Riduciamo il raggio della regione di confidenza
        options.rhoend = eps;

%       Eseguiamo di nuovo il problema di soddisfacimento partendo però da
%       un ansatz già buono
        [x_opt, ~, ~, output] = prima([], x_init, [], [], [], [], lb, ub, nonlcon, options);

%       Aggiorniamo il contatore delle iterazioni
        it_opt = it_opt + 1;
    end

%   Controlliamo se abbiamo raggiunto o meno la precisione desiderata
    if (output.constrviolation < 1e-5)
        pass = true;
    end

%   Recuperiamo le codewords
    x_tot = zeros(2*dim,1);

    x_tot(supp0) = x_opt(1:dim/2);
    x_tot(supp1+dim) = x_opt(dim/2+1:end);

    cw(:,1) = x_tot(1:dim);
    cw(:,2) = x_tot(dim+1:end);

%   Genriamo anche le errorwords
    [ew0, ew1] = genera_errorwords(cw);

%   Output per sapere l'esito della ricerca
    output.constrviolation
end