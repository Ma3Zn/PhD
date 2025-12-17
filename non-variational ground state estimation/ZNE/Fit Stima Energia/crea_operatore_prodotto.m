function op = crea_operatore_prodotto(spin_chain, dim_vals, pos, id)
%   Funzione che crea l'opportuno operatore di spin per la posizione e l'id
%   fornito

%   Dimensione dello spazio
    dim = 1;

%   Calcoliamo la dimensione dello spazio
    for it = 1:length(spin_chain)
        dim = dim * dim_vals(it);
    end

%   Inizializziamo l'operatore nel modo opportuno
    if(pos == 1)
        op = sparse(sop(spin_chain(1), id));
    else
        op = speye(dim_vals(1));
    end

%   Cicliamo sugli spin del sistema
    for it = 2:length(spin_chain)
        if (pos == it)
            op = kron(op, sparse(sop(spin_chain(pos), id)));
        else
            op = kron(op, speye(dim_vals(it)));
        end
    end
end