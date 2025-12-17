function lambda = calcola_lambda(j, in)
%   Funzione che dati in input la rho attuale ed i parametri del sistema
%   calcola i rispettivi lambda e li ritorna in un vettore riga

%   Allochiamo lambda
    lambda = zeros(1,in.dim^2);

%   Cicliamo sulle varie rho_k su cui calcolare la traccia
    for n = 1:in.dim
        for m = 1:in.dim
%           Recuperiamo l'indice iterativo attuale
            k = (n-1) * in.dim + m;

%           Calcoliamo l'evoluzione del sistema quantistico dovuta alla
%           nostra quantum operatrion
%
%           Data la forma peculiare di rho_j possiamo evitare di eseguire
%           il prodotto matrice vettore, in quanto ...
            rho = devettorizza_matrice(in.evol(:,j));

%           Calcoliamo il valore di lambda_jk come la traccia tra la rho
%           ottenuta dall'evoluzione dovuta alla quantum operation e rho_k'
%           Poichè rho_k' = |m><n| --> la traccia desiderata è = rho(n,m)
            lambda(k) = rho(n,m);
        end
    end
end