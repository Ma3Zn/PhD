function CB_completa = crea_base_completa(CB, q, dim_a)
%   Funzione che dati in input la matrice del cambiamento di base e
%   l'indice del qudit su cui applicare un'effettiva EC calcola la matrice
%   del cambiamento di base per il sistema completo

    switch q
        case 1
            CB_completa = kron(CB, speye(dim_a));
            CB_completa = kron(CB_completa, kron(CB, CB));
        case 2
            CB_completa = kron(CB, CB);
            CB_completa = kron(CB_completa, kron(speye(dim_a), CB));
        case 3
            CB_completa = kron(CB, CB);
            CB_completa = kron(CB_completa, kron(CB, speye(dim_a)));
        otherwise
            disp('ERRORE:: indice del qudit da correggere sbagliato');
    end
end