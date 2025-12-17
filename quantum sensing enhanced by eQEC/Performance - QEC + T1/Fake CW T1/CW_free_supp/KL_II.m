function [KL] = KL_II(x)
%   Funzione che calcola la seconda condizione di Knill-Laflamme.
%   x è un vettore di dimensione 2*d contenente gli elementi della codifica
%   di |0_L> nelle prime d entrate e quelli della codifica di |1_L> nelle
%   restanti d.
%   Necessità di una struttura globale ERR contenete gli errori del sistema
%   necessaria per non dover codificare a mano gli errori all'interno della
%   funzione

    global ERR;

    n_err = ERR.n;
    dim   = ERR.dim;

    KL = zeros(n_err^2, 1);

    for j = 1:n_err
        for k = 1:n_err
%           Calcoliamo l'indice attuale
            idx = (j-1) * n_err + k;

%           Calcoliamo il valore delle condizioni di KL
            KL(idx) = x(1:dim)' * ERR.E{j}' * ERR.E{k} * x(dim+1:2*dim);

%           Trasformiamo le condizioni puramente immaginarie in condizioni
%           reali
            if (~isreal(KL(idx)) && real(KL(idx)) < (1e2 * eps) )
                KL(idx) = imag(KL(idx))
            end
        end
    end
end