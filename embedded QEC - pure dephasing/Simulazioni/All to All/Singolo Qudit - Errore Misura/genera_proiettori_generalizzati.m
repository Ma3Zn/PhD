function [prj_gen, p] = genera_proiettori_generalizzati(prj, dim, p_err_mis)
%   Creiamo i proiettori generalizzati ("storti") del codice, per includere
%   anche un errore di misura nelle simulazioni

    prj_gen = cell(1,dim/2);
    p = cell(1,dim/2);

    for k = 1:dim/2
%       Generiamo la distribuzione di errore sui vari k, i.e., la
%       "stortezza" dei vari proiettori
        p{k} = genera_distribuzione_errore_misura(p_err_mis, dim/2, k);

        prj_gen{k} = sparse(dim^2/2, dim^2/2);

        for j = 1:dim/2
            prj_gen{k} = prj_gen{k} + p{k}(j) * prj{j};
        end
    end
end