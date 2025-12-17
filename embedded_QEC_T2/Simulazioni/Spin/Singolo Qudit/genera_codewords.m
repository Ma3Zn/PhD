function [M] = genera_codewords(spin)
%   Funzione che dato uno spin genera le codewords binomiali

%   Dimensione dello spazio
    dim = 2*spin + 1;

%   Generiamo l'operatore d'errore
    Sz = diag(-spin:spin);

%   Matrice che conterrà in ogni colonna le codewords relative, raggruppate
%   a blocchi in base ad l (=0,1)
    M = zeros(dim,dim);

    [M(:,1), M(:,dim/2 + 1)] = codewords_binomiali(spin);

    for j = 2:dim/2
        tmp = Sz * M(:,j - 1);
        tmp = tmp / sqrt(tmp'*tmp);

        M(:,j) = gs(M(:,1 : j - 1), tmp);

        tmp = Sz * M(:,j - 1 + dim/2);
        tmp = tmp / sqrt(tmp'*tmp);

        M(:,j + dim/2) = gs(M(:,1 + dim/2 : j - 1 + dim/2), tmp);
    end
end