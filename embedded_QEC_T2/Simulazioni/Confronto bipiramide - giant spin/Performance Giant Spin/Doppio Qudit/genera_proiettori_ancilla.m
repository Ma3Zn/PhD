function prj = genera_proiettori_ancilla(dim)
%   Fuzione che dato uno spin genera i proiettori per il codice binomiale
%   agli stabilizzatori

    prj = cell(1,dim/2);
    tmp = speye(dim);

    for j = 1:dim/2
        tmp_1 = sparse(dim/2,dim/2);
        tmp_1(j,j) = 1;

        prj{j} = kron(tmp, tmp_1);
    end
end