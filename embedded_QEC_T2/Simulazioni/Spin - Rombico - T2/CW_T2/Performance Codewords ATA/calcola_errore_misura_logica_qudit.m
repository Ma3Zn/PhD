function err = calcola_errore_misura_logica_qudit(sim, rho)
%   Funzione che calcola l'errore come l'errore che compiremmo eseguendo una 
%   misura logica del qudit al termine di un effettivo circuito

%   Recuperiamo la rho finale
    rho_fin = sim.psi_fin * sim.psi_fin';

%   Calcoliamo le probabilità di misura esatte dei vari stati logici
    p_esatta(1) = real(trace(rho_fin * sim.prj_l{1}));
    p_esatta(2) = real(trace(rho_fin * sim.prj_l{2}));

%   Calcoliamo le probabilità di misura esatte dei vari stati logici
    p_inc(1) = real(trace(rho * sim.prj_l{1}));
    p_inc(2) = real(trace(rho * sim.prj_l{2}));

    err = abs(p_esatta(2) - p_esatta(1) - (p_inc(2) - p_inc(1)));
end