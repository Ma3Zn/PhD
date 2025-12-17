function ff = aumenta_ff(folding_factors, N)
%   Funzione che per ogni valore di folding_factors genera un vettore di N
%   valori identici

%   Recuperiamo la lunghezza di folding_factors
    l = length(folding_factors);

%   Inizializziamo il vettore contenente gli ff finali
    ff = zeros(1, l * N);

%   Generiamo un vettore di 0 della dimensione opportuna
    z = zeros(1, N);

%   Cicliamo sui folding factors ed ampliamo il vettore nel modo opportuno
    for i = 1:l
        ff((i-1)*N + 1 : i*N) = z + folding_factors(i);
    end
end 