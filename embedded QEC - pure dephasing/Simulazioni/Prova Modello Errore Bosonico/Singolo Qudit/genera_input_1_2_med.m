function in = genera_input_1_2_med(T2, ang, amp, T)
%   Funzione che dati in input opportuni parametri genera l'input per le
%   simulazioni delle prestazioni del codice di EC mediate su diverse
%   operazioni logiche e diversi stati iniziali

%   Creiamo il gate fisico da simulare
    gate_fisico = rotazione_piana(ang(1), ang(2));

% % % %   Memory Time
% % %     gate_fisico = eye(2);

%   Calcoliamo il logaritmo dell'unitaria che vogliamo implementare
    log_op = logm(gate_fisico);

%   Generiamo la parte coerente
    log_sop_evol = genera_superoperatore_sx(log_op) ...
                 - genera_superoperatore_dx(log_op);

%   Operatore d'errore (opportunamente ordinato)
    E = [0 1; 0 0];

%   Creiamo la parte incoerente
    rateo = 1/T2 * T;

    sop_E = 2*genera_superoperatore_sx(E) * genera_superoperatore_dx(E')...
            - genera_superoperatore_sx(E'*E) ...
            - genera_superoperatore_dx(E'*E);
    sop_E = rateo * sop_E;

%   Aggiungiamo la parte incoerente
    log_sop_evol = log_sop_evol + sop_E;

%   Esponenziamo l'operaore
    sop_evol = expm(log_sop_evol);

%   rho iniziale della simulazione
    psi_in = amp(1) * [1 0]' + amp(2) * [0 1]';
    rho_in = psi_in*psi_in';

%   Vettorizziamo rho
    rho = vettorizza_matrice(rho_in);

%   Stato finale del sistema
    psi_fin = gate_fisico * psi_in;

    in.sop_evol = sop_evol;
    in.psi_fin = psi_fin;
    in.rho = rho;
end