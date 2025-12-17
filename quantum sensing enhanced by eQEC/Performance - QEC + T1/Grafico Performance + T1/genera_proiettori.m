function [prj] = genera_proiettori(CB, dim)
%   Fuzione che data la dimensione del qudit genera i proiettori di
%   stabilizzazione dell'ancilla

    prj = cell(1,dim/2);
    tmp = speye(dim);

    for j = 1:dim/2
        tmp_1 = sparse(dim/2,1);
        tmp_1(j) = 1;
        tmp_1 = tmp_1*tmp_1';

        prj{j} = kron(tmp, tmp_1);
    end
end