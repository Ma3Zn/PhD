function rho_rid = traccia_ancille(sim, rho, q)
%   Funzione che dato un sistema completo esegue tutte le operazioni
%   necessarie per tracciare fuori tutte le ancille.

    I = speye(sim.dim_q^2);

%   Recuperiamo la matrice di riordinamento opportuna
    P = sim.P;

%   Eseguiamo le opportune operazioni in base a dove sta l'ancilla
    switch q
        case 1
%           Calcoliamo l'opportuno riordinamento della base
            riord = kron(P, I);
            
%           Riordiniamo la base scambiando ancilla e control [QAQS] --> [AQQS]
            rho = riord' * rho * riord;

%           Tracciamo fuori l'ancilla relativa al control [AQQS] --> [QQS]
            rho_rid = partial_trace_0(rho, sim.dim_a, sim.dim_q^3);
        case 2
%           Calcoliamo l'opportuno riordinamento della base
            riord = kron(I, P);

%           Riordiniamo la base scambiando ancilla e control [QQAS] --> [QQSA]
            rho = riord * rho * riord';

%           Tracciamo fuori l'ancilla relativa al control [QQSA] --> [QQS]
            rho_rid = partial_trace_1(rho, sim.dim_q^3, sim.dim_a);
        case 3
%           Tracciamo fuori l'ancilla relativa allo switch [QQSA] --> [QQS]
            rho_rid = partial_trace_1(rho, sim.dim_q^3, sim.dim_a);
        otherwise
            disp('ERRORE:: indice del qudit da correggere sbagliato');
    end
end