function E = incoherent_krauss_operators(in, Err)
%   Funzione che dati gli operatori di krauss derivanti da una tomografia
%   di un processo, ricava gli operatori di krauss che caratterizzano la
%   sola parte incoerente. Restituisce tutti e soli gli operatori che hanno
%   norma (di frobenius) maggiore di tol.

%   Recuperiamo il numero di operatori di Krauss
    sz = size(Err);

%   Allochiamo gli operatori di krauss incoerenti
    E = zeros(in.dim, in.dim, sz(3));

%   Rimuoviamo la parte coerente
    for j = 1:sz(3)
        E(:,:,j) = in.U'*Err(:,:,j);
    end

end