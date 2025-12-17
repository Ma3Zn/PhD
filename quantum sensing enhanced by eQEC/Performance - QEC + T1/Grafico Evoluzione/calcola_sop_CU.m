function sop_ev = calcola_sop_CU(sim, perfect_EC)
%   Funzione che genera il superoperatore di evoluzione incoerente del CU
%   di durata sim.durata_op_EC

%   Matrice identità di suporto
    I = speye(sim.dim_a);

%   costruiamo la forma base del CU
    CU = genera_CU(sim.dim_q);

%   Calcoliamo il logaritmo di questo operatore
    lCU = logm(CU);

%   scriviamo il logaritmo in base computazionale
    lCU = kron(sim.CB,I) * lCU * kron(sim.CB,I)';
        
%   Rimuoviamo eventuale sporcizia numerica
%     lCU = pulisci_matrice(lCU, 1e3 * eps);

%   Calcoliamo il superoperatore commutatore
    sop_comm_CU = genera_superoperatore_sx(lCU) ...
                - genera_superoperatore_dx(lCU) ;

%   Aggiungiamo la parte incoerente
    sop_comm_CU = sop_comm_CU + sim.sop_inc_completo * (sim.durata_op_EC / sim.T2);

%   Calcoliamo il superoperatore di evoluzione
    sop_ev = fastExpm(sop_comm_CU);
end