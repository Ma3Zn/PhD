function sop_evol = genera_operatore_evoluzione_QO(in)
%   Funzione che dato in input i parametri necessari genera il
%   superoperatore d'evoluzione caratterizzante la nostra quantum operation

%   Calcoliamo il logaritmo dell'unitaria che vogliamo implementare
    op = sparse(logm(in.U));

%   Generiamo la parte coerente
    sop_evol = genera_superoperatore_sx(op) ...
             - genera_superoperatore_dx(op);

% %   Aggiungiamo la parte incoerente
%     sop_evol = sop_evol + (in.t / in.T1) * in.sop_inc_T1 ...
%                         + (in.t / in.T2) * in.sop_inc_T2;

% %   Solo T1
%     sop_evol = sop_evol + (in.t / in.T1) * in.sop_inc_T1;

%   Solo T2
    sop_evol = sop_evol + (in.t / in.T2) * in.sop_inc_T2;

%   Esponenziamo l'operaore
    sop_evol = fastExpm(sop_evol);
end