function [E, pesi] = normalizza_operatori(ops)
%   Funzione che dati in input un array di operatori li normalizza secondo
%   la norma di frobenius

%   Recuperiamo le dimensioni degli operatori
    sz = size(ops);

%   Allochiamo gli operatori normalizzati
    E = zeros(sz);

%   Allochiamo i pesi
    pesi = zeros(sz(3),1);

%   Normalizziamo gli operatori
    for i = 1:sz(3)
        pesi(i)  = norm(ops(:,:,i), 'fro');

%       Verifichiamo che l'operatore non sia un operatore nullo
        if pesi(i) > 0
            E(:,:,i) = 1/pesi(i) * ops(:,:,i);
        end
    end
end