function k = seleziona_indice(vot)
%   Funzione che date le votazioni per gli indici di recovery ritorna
%   l'esito della votazione per maggioranza dell'indice più votato. Se non
%   si raggiunge una maggioranza assoluta si ritorna l'ultimo indice più
%   votato

%   Recuperiamo il numero di votazioni
    n_vot = length(vot);

%   Recuperiamo gli indici votati
    idxs = unique(vot);

%   Recuperiamo il numero di ripetizioni per ogni indice votato
    rep_vot = histc(vot, idxs);

%   Recuperiamo il massimo numero di votazioni
    max_rep_vot = max(rep_vot);

%   Recuperiamo gli indici degli indici con conteggio maggiore
    i_max = find(rep_vot == max_rep_vot);

%   Calcoliamo quanti indici sono stati votati 'max_n_vot' volte
    n_idx_max = histc(rep_vot, max_rep_vot);

    if (n_idx_max == 1)

        k = idxs(i_max);

    else

        k = min(idxs(i_max));

%         j = n_vot;
%         while j > 0
% 
%             jj = n_idx_max;
%             while jj > 0
% 
%                 if (vot(j) == idxs(i_max(jj)))
%                     
%                     k = idxs(jj);
% 
%                     j = 0;
%                     jj = 0;
%                 end
% 
%                 jj = jj - 1;
%             end
% 
%             j = j - 1;
%         end

    end
end