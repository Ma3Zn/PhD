function prj = genera_proiettori_ideali(dim)
%   Fuzione che dato uno spin genera i proiettori per il codice binomiale
%   agli stabilizzatori

    prj = cell(dim/2);

    for j = 1:dim/2
        tmp = sparse(dim,dim);
        tmp(j,j) = 1;
        tmp(j + dim/2, j + dim/2) = 1;

        prj{j} = tmp;
    end
end