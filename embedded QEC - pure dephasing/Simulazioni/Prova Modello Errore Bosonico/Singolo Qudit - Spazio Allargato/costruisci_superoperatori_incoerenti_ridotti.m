function sop_err = costruisci_superoperatori_incoerenti_ridotti(sim)
%   Funzione che genera i superopertori incoerenti del sistema dati in
%   input gli opportuni parametri

%   Recuperiamo l'operatore d'errore del sistema
    op_err = sparse(diag(sqrt(1:sim.dim_q-1),1));

%   Generiamo il superoperatore del lindbladiano per l'evoluzione
%   incoerente del sistema (attenzione qui usiamo il fatto che Sz sia
%   hermitiano)
    sop_err = 2*genera_superoperatore_sx(op_err) * genera_superoperatore_dx(op_err')...
            - genera_superoperatore_sx(op_err' * op_err) ...
            - genera_superoperatore_dx(op_err' * op_err);
end