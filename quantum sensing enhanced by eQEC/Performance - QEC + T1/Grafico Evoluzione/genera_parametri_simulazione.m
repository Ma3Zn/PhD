function sim = genera_parametri_simulazione(S, Bx, T2, T_rot)
%   Funzione che genera i parametri per la simulazione Monte Carlo del
%   protocollo di sensing

%   Generiamo il sistema fisico di base
    sim = genera_sistema(S, Bx, T2);

%   Generiamo i superoperatori incoerenti del sistema
    sim.sop_inc_completo = calcola_sop_inc_completo(sim);

%   Generiamo i superoperatori incoerenti ridotti del sistema
    sim.sop_inc_ridotto = calcola_sop_inc_ridotto(sim);

%   Stabiliamo la durata delle operazioni di stabilizzazione e recovery
    sim.durata_op_EC = sim.t_es(sim.dim_q);

%   Stabiliamo la durata della misura (review + 50ns detection fotone)
    sim.durata_misura = 300;

%   Generiamo il superoperatore di evoluzione per il CU
    sim.sop_CU = calcola_sop_CU(sim);

%   Generiamo i superoperatori di evoluzione per la recovery
    sim.sop_R = calcola_sop_R(sim);

%   Generiamo il superoperatore di idle per simulare la lunghezza del
%   processo di misura
    sim.sop_idle = calcola_sop_idle(sim);

%   Generiamo il superoperatore del commutatore del generatore della
%   rotazione logica
    [sim, sim.sop_comm_rot] = calcola_sop_comm_rot(sim);

%   Settiamo per quanto tempo far evolvere il sistema con una rotazione
%   logica tramite due applicazioni successive di error correction
    sim.dt = T_rot;

%   Durata intera procedura EC
    sim.durata_EC = 2 * sim.durata_op_EC + sim.durata_misura;

%   Durata totale di uno step di rotaione logica e di EC
    sim.dt_tot = sim.dt + sim.durata_EC;

%   Generiamo il superoperatore di evoluzione per un pezzo della rotazione
%   logica
    sim.sop_rot = calcola_sop_rot(sim);

%   Tempo di evoluzione sotto rotazione logica
    sim.T = 0;

%   Contatore cicli QEC eseguiti
    sim.it = 0;

%   Tempo finale di simulazione
    sim.n = 1e5;
    sim.T_fin = sim.dt_tot * sim.n; % [ ns ]

%   Allochiamo lo spazio per il salvataggio dell'evoluzione della
%   popolazione
    sim.pop = zeros(1,sim.n+1);
end
