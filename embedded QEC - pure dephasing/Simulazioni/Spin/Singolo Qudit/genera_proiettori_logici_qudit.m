function prj_l = genera_proiettori_logici_qudit(dim, dim_q)
%   Funzione che crea i proiettori logici per la misura logica del qubit
%   virtualizzato

    prj_l = cell(1,2);

    I = eye(dim_q/2);
    Z = sparse(dim_q/2, dim_q/2);

    prj_l{1} = [I Z; Z Z];
    prj_l{2} = [Z Z; Z I];

end