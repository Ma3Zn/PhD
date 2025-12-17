function sop_ev = calcola_sop_idle(sim)
%   Funzione per il calcolo del superoperatore di idle sul sistema completo
%   (qudit + ancilla) per simulare la lunghezza del processo di misura

%   scriviamo il logaritmo in base computazionale
    l_idle = sparse(sim.dim, sim.dim);

%   Calcoliamo il superoperatore commutatore
    sop_comm_idle = genera_superoperatore_sx(l_idle) ...
                  - genera_superoperatore_dx(l_idle) ;

%   Aggiungiamo la parte incoerente
    sop_comm_idle = sop_comm_idle + sim.sop_inc_completo * (sim.durata_misura / sim.T2);

%   Calcoliamo il superoperatore di evoluzione
    sop_ev = fastExpm(sop_comm_idle);
end