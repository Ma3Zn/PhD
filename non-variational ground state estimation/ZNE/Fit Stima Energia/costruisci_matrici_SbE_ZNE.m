function [H, S] = costruisci_matrici_SbE_ZNE(exp_val)
% Funzoine che passati i valori di aspettazione di N potenze
% dell'Hamiltoniana ne costruisce le relative matrici di una Krylov-like
% SbE

%   Aggiungiamo anche il valore di aspettazione dell'identità
    exp_val = [1 exp_val];

%   Recuperiamo il numero di potenze dell'Hamiltoniana da usare
    N = length(exp_val) - 1;

%   Calcoliamo la dimensione della SbE [il +1 è perchè usiamo anche I] 
    dim = (N + 1)/2;

%   Allochiamo le matrici della SbE
    H = zeros(dim);
    S = zeros(dim);

%   Costruiamo le matrice della SbE
    idx = 1;
    l_vec = 0;
    for n_diag = -(N-1)/2 : (N-1)/2

%       Calcoliamo la lunghezza opportuna del vettore di elementi
        if n_diag <= 0
            l_vec = l_vec + 1;
        else
            l_vec = l_vec - 1;
        end

%       Costruiamo il vettore di elementi
        vec = ones(1,l_vec);

%       Aggiorniamo H
        H = H + diag(exp_val(idx+1) * vec, n_diag);

%       Aggiorniamo S
        S = S + diag(exp_val(idx) * vec, n_diag);

%       Aggiorniamo l'indice di avanzamento
        idx = idx + 1;
    end

%   Ordianiamo in modo opportuno le due matrici della SbE
    H = H(end:-1:1,:);
    S = S(end:-1:1,:);
end