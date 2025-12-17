function [errore_medio, T] = calcola_errore_medio(dim,T2_vals)
%   Funzione che dati in input uno spin S ed i valori di T2 da testare,
%   calcola la fidelity della procedura di EC agli stabilizzatori eseguendo
%   una completa decomposizione della procedure richieste.
%
%   Fornisce in output l'errore medio commesso dalla procedura al variare
%   di T2

%   Numero di T2 da provare
    N_T2 = length(T2_vals);

%   Allocazione variabile di output
    errore_medio = zeros(N_T2,1);
    T            = zeros(N_T2,1);

%   Teniamo traccia dell'indice di iterazione
    idx = 1;

%   Generazione input necessario all'esecuzione del circuito
    sim = genera_input_QEC(dim);

%   Recuperiamo la rho iniziale
    rho_in = sim.rho;

    for T2 = T2_vals

%       Resettiamo la rho iniziale [inutle perchè questa rho non viene mai modificata]
        sim.rho = rho_in;

%       Op evoluzione incoerente relativo al gate logico
        log_G  = logm(sim.gate_logico);
        log_G  = pulisci_matrice(log_G, 1e3 * eps);
        comm_G = genera_superoperatore_sx(log_G) - genera_superoperatore_dx(log_G);

%       Calcoliamo il superoperatore di evoluzioine
        comm_G = comm_G + sim.sop_inc_completo * (30/T2);
        sim.ev_gl = fastExpm(comm_G);

%       OP evoluzione incoerente relativo al CU
        comm_CU = sim.sop_comm_CU;
        comm_CU = comm_CU + sim.sop_inc_completo * (30/T2);
        sim.ev_CU = fastExpm(comm_CU);

%       OP evoluzione incoerente relativo alle operazioni di recovery
        sim.ev_R = genera_operatore_evoluzione_R(sim, 30, T2);

%       Calcoliamo l'errore del circuito implementante l'algoritmo di QEC
        [p_osservazione, errore, t] = calcola_errore(sim);

%       Calcoliamo l'errore medio della procedura
        errore_medio(idx) = p_osservazione.' * errore;
        T(idx) = t;

%       Aggiorniamo l'indice di iterazione
        idx = idx + 1;
    end
end