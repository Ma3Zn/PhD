function in = genera_input_1_2(T2, ang, amp_c, amp_t, T)
%   Funzione che dati in input opportuni parametri genera l'input per le
%   simulazioni delle prestazioni del codice di EC mediate su diverse
%   operazioni logiche e diversi stati iniziali

%   Creiamo il gate fisico da simulare
    gate_fisico = [eye(2) zeros(2); zeros(2) rotazione_Z(ang)];
    gate_fisico = pulisci_matrice(gate_fisico, 1e2 * eps);

%   Calcoliamo il logaritmo dell'unitaria che vogliamo implementare
    log_op = logm(gate_fisico);

%   Generiamo la parte coerente
    log_sop_evol = genera_superoperatore_sx(log_op) ...
                 - genera_superoperatore_dx(log_op);

%   Operatore d'errore (opportunamente ordinato)
    Ec = kron(1/2 * [-1 0; 0 1], eye(2));
    Et = kron(eye(2), 1/2 * [-1 0; 0 1]);

%   Creiamo la parte incoerente
    rateo = 1/T2 * T;

    sop_Ec = 2*genera_superoperatore_sx(Ec) * genera_superoperatore_dx(Ec)...
            - genera_superoperatore_sx(Ec*Ec) ...
            - genera_superoperatore_dx(Ec*Ec);
    sop_Ec = rateo * sop_Ec;

    sop_Et = 2*genera_superoperatore_sx(Et) * genera_superoperatore_dx(Et)...
            - genera_superoperatore_sx(Et*Et) ...
            - genera_superoperatore_dx(Et*Et);
    sop_Et = rateo * sop_Et;

%   Aggiungiamo la parte incoerente
    log_sop_evol = log_sop_evol + sop_Ec + sop_Et;

%   Esponenziamo l'operaore
    sop_evol = expm(log_sop_evol);

%   rho iniziale della simulazione
    psi_in_c = amp_c(1) * [1 0]' + amp_c(2) * [0 1]';
    psi_in_t = amp_t(1) * [1 0]' + amp_t(2) * [0 1]';
    psi_in = kron(psi_in_c, psi_in_t);
    rho_in = psi_in*psi_in';

%   Vettorizziamo rho
    rho = vettorizza_matrice(rho_in);

%   Stato finale del sistema
    psi_fin = gate_fisico * psi_in;

    in.sop_evol = sop_evol;
    in.psi_fin = psi_fin;
    in.rho = rho;
end