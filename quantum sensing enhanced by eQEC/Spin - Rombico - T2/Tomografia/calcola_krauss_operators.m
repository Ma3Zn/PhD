function E = calcola_krauss_operators(U, d, in)
%   Funzione che dati gli autovettori di CHI ed i rispettivi autovalori
%   calcola la rappresentazione degli operatori Krauss della quantum
%   operation

%   Allochiamo E
    E = cell(1,in.dim^2);

%   Cicliamo sull'indice dell'operatore di Krauss
    for i = 1:in.dim^2
        E{i} = sparse(in.dim, in.dim);

        for j = 1:in.dim^2
            E{i} = E{i} + U(j,i) * in.A{j};
        end

        E{i} = pulisci_matrice(E{i} * sqrt(d(i)), 1e4 * eps);
    end
end