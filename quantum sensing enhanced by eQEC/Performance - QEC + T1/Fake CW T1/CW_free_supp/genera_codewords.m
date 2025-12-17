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
    [E, n] = operatori_errori(dim);

%   Generazione parametri errore
    global ERR;

%   Calcoliamo le dimensioni del sistema
    ERR.dim  = dim;

%   Recuperiamo TUTTI gli operatori d'errore rilevanti
    ERR.E = cell(1,ERR.n);

%   Inseriamo gli operatori di Krauss che vogliamo correggere
    for j = 1:ERR.n
        ERR.E{j} = E{j};
    end

%   Inseriamo lo spacing per le codewords. Sempre 1 per T2
    ERR.spacing = 3;

%   Creazione handle per la funzione da ottimizzare e per i vincoli non
%   lineari
    nonlcon = @nonlinear_constraint;

%   Generiamo upper e lower bound per i valori delle variabili di encoding
    lb = -1 * ones(2*ERR.dim,1);
    ub =      ones(2*ERR.dim,1);

%   Generiamo la condizione iniziale
%     if dim < 13
%         x_init = genera_condizione_iniziale_STD(dim);
%     else
        x_init = genera_condizione_iniziale_MOD(dim);
%     end

%   Settiamo le opzioni di ricerca per il solver
    options = struct();
    options.maxfun = 5e2*length(x_init);

%   Eseguiamo il primo step di ottimizzazione
    [x_opt, ~, ~, output] = prima([], x_init, [], [], [], [], lb, ub, nonlcon, options);

%   Contatore per il numero di iterazioni del ciclo di ottimizzazione
    it_opt = 1;

% % % %   Iteriamo l'ottimizzazine
% % %     while (output.constrviolation > tol && it_opt < 10)
% % % %       Recuperiamo le codewords
% % %         cw(:,1) = x_opt(1:dim);
% % %         cw(:,2) = x_opt(dim+1:end);
% % % 
% % % %       forziamo i vincoli di spacing sulle codewords [NO per T1]
% % %         cw = pulisci_codewords(cw);
% % % 
% % % %       Aggiorniamo la condizione iniziale con le codewords trovate al
% % % %       passo precedente
% % %         x_init = [cw(:,1); cw(:,2)];
% % % 
% % % %       Riduciamo il raggio della regione di confidenza
% % %         options.rhoend = eps;
% % % 
% % % %       Eseguiamo di nuovo il problema di soddisfacimento partendo però da
% % % %       un ansatz già buono
% % %         [x_opt, ~, ~, output] = prima([], x_init, [], [], [], [], lb, ub, nonlcon, options);
% % % 
% % % %       Aggiorniamo il contatore delle iterazioni
% % %         it_opt = it_opt + 1;
% % %     end

%   Controlliamo se abbiamo raggiunto o meno la precisione desiderata
    if (output.constrviolation < 1e-5)
        pass = true;
    end

%   Recuperiamo le codewords
    cw(:,1) = x_opt(1:dim);
    cw(:,2) = x_opt(dim+1:end);

% % % %   Forziamo i vincoli di spacing sulle codewords
% % %     cw = pulisci_codewords(cw);

%   Genriamo anche le errorwords
    [ew0, ew1] = genera_errorwords(cw);

%   Output per sapere l'esito della ricerca
    output.constrviolation
end