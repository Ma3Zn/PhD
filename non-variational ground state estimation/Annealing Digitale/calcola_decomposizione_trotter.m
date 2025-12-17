function U = calcola_decomposizione_trotter(trt, Ha, Hb, dt)
%   Funzione per il calcolo dell'operatore di evoluzione dato dalla
%   trotterizzazione opportuna

%   Recuperiamo la dimensione del sistema
    dim = length(Ha);

%   Recuperiamo il numero di cicli da eseguire
    n = length(trt.a);

%   Allochiamo lo spazio per le varibili necessarie
    U  = speye(dim);

%   Cicliamo sul numero di coefficienti ed applichiamo la formula di
%   trotter
    for j = 1:n
        Ua = fastExpm(-1i * Ha * trt.a(j) * dt);
        Ub = fastExpm(-1i * Hb * trt.b(j) * dt);
        U = U * Ua * Ub;
    end
end