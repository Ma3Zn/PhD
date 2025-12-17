function E = riordina_operatori(E, n_op)
%   Funzione che riordina gli operatori di Krauss in base alla loro norma

%   Calcoliamo le norme degli operatori
    norma = zeros(n_op,1);

    for i = 1:n_op
        norma(i) = norm(full(E{i}), 'fro');
    end

%   ordianiamo le norme degli operatori
    [~,I] = sort(norma, 'descend');

%   Riordiniamo gli operatori
    E = E(I);

end