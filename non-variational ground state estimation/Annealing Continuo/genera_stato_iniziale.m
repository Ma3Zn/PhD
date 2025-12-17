function psi = genera_stato_iniziale(sim)
%   Funzione che dati i parametri della simulazione genera lo stato
%   iniziale di questa

%   Generiamo l'Hamiltoniana iniziale
    H0 = sim.H(sim, 0);

%   Calcoliamo gli autovettori dell'Hamiltoniana
    [V, ~] = eig(full(H0));

%   Inizializziamo il sistema nel ground state iniziale
    psi = V(:,1);
end