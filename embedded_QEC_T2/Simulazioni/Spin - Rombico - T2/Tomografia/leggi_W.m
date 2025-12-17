function W = leggi_W(dim)
%   Funzione per la lettura dei parametri per il T1

%   Leggiamo i ratei di T1
    W = load("par/W.txt");

%   Scaliamo opportunamente i parametri
    W = W * 1e2;

%   Rimuoviamo la diagonale di W che non è necessaria (equivale a sommare
%   su b \ne a nel lindbladiano).
    W = W - diag(diag(W));

%   Ci serve la trasposta di W
    W = W.';

%   Riduciamo le dimensioni di W
    W = W(1:dim, 1:dim);
end