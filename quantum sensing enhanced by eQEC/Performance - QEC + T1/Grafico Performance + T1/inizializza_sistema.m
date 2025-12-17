function sim = inizializza_sistema(sim)
%   Funzione per l'inizializzazione del sistema per la prossima esecuzione 
%   del protocollo 

%   Ripristiniamo il dt iniziale
    sim.dt = 2 * sim.durata_op_EC;

%   Ricostruiamo il superoperatore di evoluzione della rotazione logica
    sim.sop_rot = calcola_sop_rot(sim);

%   Ripristiniamo il tempo stimato del minimo
    sim.T = 0;

%   Ripristiniamo lo stato iniziale del sistema
    psi_in    = zeros(sim.dim_q,1);
    psi_in(1) = 1;
    
    psi_in = sim.CB * psi_in;

    sim.rho = vettorizza_matrice(psi_in * psi_in');

    sim.rho_old = sim.rho;
end