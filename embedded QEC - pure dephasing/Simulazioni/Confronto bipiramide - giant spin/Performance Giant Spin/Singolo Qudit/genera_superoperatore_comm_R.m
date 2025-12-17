function [R, sop_comm_R] = genera_superoperatore_comm_R(dim, CB)
%   Funzione che dati in input la dimensione degli oggetti del sistema e la
%   matrice del cambio di base [Logica --> Computazionale] per il sistema
%   ridotto calcola gli operatori di recovery ed i relativi superoperatori
%   di commutazione in base computazionale

%   Generiamo gli operatori di recovery
    R = genera_operatori_recovery(dim);

%   Allochiamo lo spazio opportuno per gli operatori di recovery totale
    sop_comm_R  = cell(1,(dim/2)^3);
    lR          = cell(1,(dim/2)^3);

    for k = 1 : dim/2

%       Calcoliamo il logaritmo degli operatori
        lR{k} = sparse(logm(R{k}));

%       riscriviamo gli operatori in base computazionale
        lR{k} = CB * lR{k} * CB';
        R{k}  = CB *  R{k} * CB';

%       Rimuoviamo eventuale sporcizia numerica
        lR{k} = pulisci_matrice(lR{k}, 1e3 * eps);
        R{k}  = pulisci_matrice( R{k}, 1e3 * eps);

%       Scriviamo il superoperatore associato
        sop_comm_R{k} = genera_superoperatore_sx(lR{k}) ...
                      - genera_superoperatore_dx(lR{k}) ;
    end
end