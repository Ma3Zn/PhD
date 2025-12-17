function x = pulisci_codewords(x)
%   Funzione che date in input codewords ottenute numericamente forza la
%   disgiunzione dei supporti e rinormalizza i vettori

%   Recuperiamo la lunghezza di x
    dim = length(x);

%   Indici 0
    idx0 = 1:dim/2;

%   Indici 1
    idx1 = dim/2+1:dim;
    
%   Normalizziamo la codewords
    x(idx0) = x(idx0) / norm(x(idx0),2);
    x(idx1) = x(idx1) / norm(x(idx1),2);
end