function idxs = calcola_indici_vettorizzati(supp, dim)
%   Funzoine che calcola gli indici delle popolazioni di uno stato logico
%   all'interno di una matrice densità vettorizzata

    tmp = zeros(dim);

    for i = supp
        tmp(i,i) = 1;
    end

    idxs = find(vettorizza_matrice(tmp));
end