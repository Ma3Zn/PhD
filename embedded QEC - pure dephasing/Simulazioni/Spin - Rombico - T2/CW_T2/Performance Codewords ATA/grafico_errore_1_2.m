% Generiamo gli input per la simulazione
genera_input_1_2;

errore = zeros(length(T2_vals),1);
idx = 1;

for T2 = T2_vals
    rho = vettorizza_matrice(rho_in);

%   Eseguiamo l'operazione sul qubit non corretto
    simula_gate_fisico;

%   Calcoliamo l'errore
    rho_fin = devettorizza_matrice(rho);
    errore(idx) = 1 - psi_fin' * rho_fin * psi_fin;

%   Incrementiamo l'indice di avanzamento della simulazione
    idx = idx + 1;
end

% Facciamo il grafico dell'errore osservato
errore_medio = errore;
fai_grafico(errore_medio, T2_vals, 1e2 * eps);

clearvars -except fig T_IDLE N_T2 T2_vals dim_vals