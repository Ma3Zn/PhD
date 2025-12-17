function dev = grafico_errore_medio(dim, T2_vals, angle_vals, init_superposition, EC_iter)
%   Funzione che dati in input la dimensione del sistema, i valori di T2
%   per cui eseguire le simulazioni, i parametri delle rotazioni e gli
%   stati iniziali su cui mediare, calcola l'errore medio dovuto
%   all'esecuzione dei vari gate logici

%   Fattore moltiplicativo per l'aumento della lunghezza delle operazioni
    tau = 10;
    fattore = 3;
    fattore_misura = 0 * fattore;

%   Numero di T2 da simulare
    [~, n_T2] = size(T2_vals);

%   Numero di operazioni da testare
    [~, n_op] = size(angle_vals);

%   Numero di stati di input da provare
    [~, n_stati_iniziali] = size(init_superposition);

%   Indice di avanzamento (T2)
    idx_T2 = 1;

%   Errore medio della procedura
    errore_medio = zeros(1, n_T2);

%   Errore massimo della procedura
    errore_massimo = zeros(1, n_T2);

%   Errore minimo della procedura
    errore_minimo = ones(1, n_T2);

%   Generiamo l'input per la simulazione
    sim = genera_input_QEC(dim);

%   Cicliamo sui valori di T2
    for T2 = T2_vals

%       Inizializiamo l'indice di stato
        idx_op = 1;

%       Inizializiamo l'errore di stato da mediare
        errore_op = zeros(1,n_op);

%       Generiamo il superoperatore di evoluzione incoerente del CU
        comm_CU = sim.sop_comm_CU;
        comm_CU = comm_CU + sim.sop_inc_completo * (sim.t_es(dim)/T2);
        sim.ev_CU = fastExpm(comm_CU);

%       Generiamo il superoperatore di evoluzione incoerente delle recovery
        sim.ev_R = genera_operatore_evoluzione_R(sim, sim.t_es(dim), T2);

%       Generiamo il superoperatore di evoluzione incoerente di una fase di
%       idle pre-misura per simulare gli effetti di una misura NON
%       istantanea
        sim.ev_idle = fastExpm(sim.sop_inc_completo * fattore_misura * tau / T2);

%       Cicliamo sui parametri delle operazioni da testare
        for ang = angle_vals

%           Scriviamo il gate logico in base computazionale
            [gl_completo, gl_ridotto] = genera_gate_logico(sim, ang);
      
%           Calcoliamo il superoperatore sparso del gate logico (pulisci
%           matrice torna già unamatrice sparsa)
            log_G  = logm(full(gl_completo));
            log_G  = pulisci_matrice(log_G, 1e3 * eps);
            comm_G = genera_superoperatore_sx(log_G) - genera_superoperatore_dx(log_G);

%           Calcoliamo il superoperatore di evoluzioine
            comm_G = comm_G + sim.sop_inc_completo * (sim.t_es(dim)/T2);
            sim.ev_gl = fastExpm(comm_G);

%           Inizializiamo l'indice di stato
            idx_stato = 1;

%           Inizializiamo l'errore di stato da mediare
            errore_stato = zeros(1,n_stati_iniziali);

%           Cicliamo sugli stati iniziali del control
            for amp = init_superposition

%               Calcoliamo le condizioni iniziali del sistema
                [sim.rho, psi_in] = crea_condizioni_iniziali(sim, amp);

% % % %   Eseguaimo opportune trasformazioni per fare la misura in una base
% % % %   logica arbitraria
% % %     psi_fin = HL * psi_fin;
% % %     psi_fin = RL * psi_fin;

%               Calcoliamo le condizioni finali logiche necessarie per la
%               tomografia del sistema
                sim.psi_fin_L = rotazione_piana(ang(1), ang(2)) * [amp(1); amp(2)];

%               Calcoliamo le condizioni finali
                sim.psi_fin = gl_ridotto * psi_in;

                if (nargin == 6)
                    sim.EC_iter = EC_iter;
                end

%               Calcoliamo l'errore del circuito implementante l'algoritmo di QEC
                [p_osservazione, errore] = calcola_errore(sim);

%               Calcoliamo l'errore medio
                errore_stato(idx_stato) = p_osservazione.' * errore;

%               Facciamo avanzare l'indice di stato
                idx_stato = idx_stato + 1;
            end

%           Recuperiamo il massimo ed il minimo valore per l'errore delle
%           operazioni
            errore_massimo(idx_T2) = max(abs([errore_stato errore_massimo(idx_T2)]));
            errore_minimo(idx_T2)  = min(abs([errore_stato errore_minimo(idx_T2)]));

%           Mediamo l'errore di stato
            errore_op(idx_op) = sum(errore_stato) / n_stati_iniziali^2;

%           Facciamo avanzare l'indice di operazione
            idx_op = idx_op + 1;
        end

%       Calcoliamo l'errore medio
        errore_medio(idx_T2) = sum(errore_op) / n_op;

%       Facciamo avanzare l'indice di T2
        idx_T2 = idx_T2 + 1;
    end

%   Facciamo il grafico della curva ottenuta
    fai_grafico(errore_medio.', T2_vals, 0e2 * eps);

%   Ritorniamo gli intorni tubolari della curva ottenura
%     dev = [errore_massimo; errore_minimo];

%   L'errore minimo ha poco senso (perchè non ha senso plottare una linea
%   in 0, non perchè sia sbagliato ...)
%     dev = errore_massimo;

% % % %   Per il plot dell'errore al variare di S
    dev = errore_medio;
end
