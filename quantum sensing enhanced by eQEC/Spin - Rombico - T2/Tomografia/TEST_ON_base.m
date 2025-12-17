clear
clc

%   Settiamo la dimensione del sistema
dim = 4;

W = leggi_W(dim);

%   generiamo la base da testare
A = genera_operatori_base(W, dim);

%   Allochiamo la matrice dei risultati
ris = zeros(dim, dim);

%   Cicliamo sugli elementi della base e testiamo l'ortonormalità della
%   base
for i = 1:dim
    for j = 1:dim
        ris(i,j) = trace(A{i}' * A{j});
    end
end

pulisci_matrice(ris, 1e2 * eps)