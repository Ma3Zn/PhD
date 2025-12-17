function [E, CHI] = process_tomography(in)
%   Funzione che dati in input i parametri necessari eseguira la process
%   tomography della quantum operation caratterizzata dall'operatore di
%   evoluzione contenuto in in.evol

%   Recuperiamo la dimensione del sistema
    dim = in.dim;

%   Allochiamo la matrice dei lambda
    lambda = zeros(dim^2);

%   Allochiamo la matrice dei beta
    beta = zeros(dim^4);

%   Ciclo parallelo per il calcolo di lambda e beta
%     parfor j = 1:dim^2
    for j = 1:dim^2
%       Rendiamo locali al lavoratore parallelo i parametri della
%       simulaizone
        in_loc = in;

%       Calcoliamo i lambda_jk (j fissato, li calcoliamo per tutti i k)
        lambda(j,:) = calcola_lambda(j, in_loc);

%       Calcoliamo i beta_{(jk)(nm)}
        beta = beta + calcola_beta(j, in_loc);
    end

%   Status simulazione
    disp('Fine ciclo parallelo');

%   Vettorizziamo lambda per poter risolvere il sistema lineare di cui CHI
%   è soluzione (l'indice di riga è quello che scorre più lenamente, dunque
%   è opportuno utilizzare la funzione già realizzata)
    lambda_vect = vettorizza_matrice(lambda);

%   Risolviamo il sistema lineare che ci ridà il vettore CHI (vedi Nielsen,
%   pag. 392)
    CHI_vect = pinv(beta) * lambda_vect;

%   Devettorizziamo CHI per ottenere la matrice da diagonalizzare. (Anche
%   in questo caso l'indice di riga è il più lento a scorrere quindi la 
%   trasformazione da fare è l'inverso di quanto fatto per la
%   vettorizzazione di lambda e dunque è lecito utilizzare la funzione già
%   realizzata)
    CHI = devettorizza_matrice(CHI_vect);

%   Rimuoviamo eventuale sporcizia numerica
    CHI = pulisci_matrice(CHI, 1e2 * eps);

%   Diagonalizziamo CHI
    [U, aut] = eig(CHI);

%   Rimuoviamo eventuale sporciza numerica
    U = pulisci_matrice(U, 1e2 * eps);

%   Estraiamo la diagonale degli autovalori
    d = pulisci_matrice(diag(aut), 1e4 * eps);

%   Rimuoviamo gli elementi negativi di d
    d(d<0) = 0;

%   Calcoliamo gli operatori di Krauss della quantum operation
    E = calcola_krauss_operators(U, d, in);
end