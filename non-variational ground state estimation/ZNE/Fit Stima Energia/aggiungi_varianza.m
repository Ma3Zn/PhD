function dist = aggiungi_varianza(exp_val, std, N_sample)
%   Funzione che dato un vettore di valori di aspettazione, una deviazione
%   standard per questi ed un numero di sample estrae genera un numero
%   opportuno di campioni seguendo la distribuzione data

%   Recuperiamo le dimensioni della matrice contenente i valori medi da
%   campionare
    [n_pow, l] = size(exp_val);

%   Inizializziamo il vettore finale dei valori campionati
    dist = zeros(n_pow, l * N_sample);

%   Cicliamo sui vari valori d'aspettazione da campionare
    for j = 1:n_pow
        for i = 1 : l
            dist(j, (i-1)*N_sample + 1 : i*N_sample) = ...
                normrnd(exp_val(j,i), std, [1 N_sample]);
        end
    end
end