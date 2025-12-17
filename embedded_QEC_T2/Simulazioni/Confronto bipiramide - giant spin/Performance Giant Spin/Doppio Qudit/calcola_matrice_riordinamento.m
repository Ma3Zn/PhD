function P = calcola_matrice_riordinamento(dim_0, dim_1)
%   Funzione che dato in input le dimensioni di due qudit calcola la
%   matrice del cambiamento di base per invertire l'ordine di questi due

    dim = dim_0 * dim_1;

%   Preallochiamo la matrice del cambio di base
    P = sparse(dim, dim);

    for k_0 = 1:dim_0
        for k_1 = 1:dim_1
            P((k_0 - 1) * dim_1 + k_1, (k_1 - 1) * dim_0 + k_0) = 1;
        end
    end
end