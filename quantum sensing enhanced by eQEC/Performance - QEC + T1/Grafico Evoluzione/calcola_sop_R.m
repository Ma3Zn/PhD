function evol_R = calcola_sop_R(sim, perfect_EC)
%   Funzione che genera i superoperatori di evoluzione incoerente degli
%   step di recovery di durata sim.durata_op_EC

%   Allochiamo lo spazio per evol_R
    evol_R = cell(1, sim.dim_a);

%   Generiamo gli operatori di recovery
    R = genera_operatori_recovery(sim.dim_q);

%   Scorriamo i vari k
    for k = 1:sim.dim_a

%       Calcoliamo il logaritmo degli operatori
        lR = sparse(logm(R{k}));

%       Riscriviamo gli operatori in base computazionale
        lR = sim.CB * lR * sim.CB';

%       Rimuoviamo eventuale sporcizia numerica
        lR = pulisci_matrice(lR, 1e3 * eps);

%       Scriviamo il superoperatore associato
        sop_comm_R = genera_superoperatore_sx(lR) ...
                   - genera_superoperatore_dx(lR) ;

%       Aggiungiamo la parte incoerente
        sop_comm_R = sop_comm_R + sim.sop_inc_ridotto * (sim.durata_op_EC / sim.T2);

%       Esponenziamo l'operatore
        evol_R{k} = pulisci_matrice(fastExpm(sop_comm_R), 1e3*eps);
    end
end