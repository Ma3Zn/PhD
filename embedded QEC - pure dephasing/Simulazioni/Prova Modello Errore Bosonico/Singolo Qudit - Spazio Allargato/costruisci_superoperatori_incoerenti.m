function sop_err = costruisci_superoperatori_incoerenti(sim)
%   Funzione che genera i superopertori incoerenti del sistema dati in
%   input gli opportuni parametri

%   Superoperatore errore
    sop_err = sparse(sim.dim^2, sim.dim^2);

%   Cicliamo sul numero di oggetti che compongono il sistema
    for i = 1:2
%       Recuperiamo l'operatore d'errore del sistema
        op_err = sim.E{i};

%       Generiamo il superoperatore del lindbladiano per l'evoluzione
%       incoerente del sistema (attenzione qui usiamo il fatto che Sz sia
%       hermitiano)
        tmp = 2*genera_superoperatore_sx(op_err) * genera_superoperatore_dx(op_err')...
              - genera_superoperatore_sx(op_err' * op_err) ...
              - genera_superoperatore_dx(op_err' * op_err);

%       Sommiamo i vari operatori d'errore
        sop_err = sop_err + tmp;
    end
end