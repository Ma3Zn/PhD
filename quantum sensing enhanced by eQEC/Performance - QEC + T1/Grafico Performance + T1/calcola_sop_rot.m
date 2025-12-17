function sop = calcola_sop_rot(sim)
%   Funzione che calcola il superoperatore di evoluzione incoerente di una
%   rotazione logica di durata sim.dt

%   Aggiungiamo la parte incoerene
    sop = sim.sop_comm_rot * sim.dt ... 
        + sim.sop_inc_T1_ridotto * (sim.dt / sim.T1) ...
        + sim.sop_inc_ridotto * (sim.dt / sim.T2);

% % % %   Solo T2
% % %     sop = sim.sop_comm_rot * sim.dt ...
% % %         + sim.sop_inc_ridotto * (sim.dt / sim.T2);

%   Calcoliamo il superoperatore di evoluzione
    sop = pulisci_matrice(fastExpm(sop), 1e3*eps);
end