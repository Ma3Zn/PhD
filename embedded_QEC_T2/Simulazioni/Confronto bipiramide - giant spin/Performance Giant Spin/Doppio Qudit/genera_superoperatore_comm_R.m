function [R_prod, sop_comm_R] = genera_superoperatore_comm_R(dim, CB)
%   Funzione che dati in input la dimensione degli oggetti del sistema e la
%   matrice del cambio di base [Logica --> Computazionale] per il sistema
%   ridotto calcola gli operatori di recovery ed i relativi superoperatori
%   di commutazione in base computazionale

%   Generiamo gli operatori di recovery
    R = genera_operatori_recovery(dim);

%   Allochiamo lo spazio opportuno per gli operatori di recovery totale
    sop_comm_R  = cell(1,(dim/2)^3);
    lR_prod     = cell(1,(dim/2)^3);
    R_prod      = cell(1,(dim/2)^3);

    for kc = 1 : dim/2
        for kt = 1 : dim/2
            for ks = 1 : dim/2
                k = (kc -1) * dim^2/4 + (kt - 1) * dim/2 + ks;

%               Costruiamo il prodotto dei vari operatori
                R_prod{k} = kron(R{kc}, kron(R{kt}, R{ks}));

%               Calcoliamo il logaritmo degli operatori
                lR_prod{k} = sparse(logm(R_prod{k}));

%               riscriviamo gli operatori in base computazionale
                lR_prod{k} = CB * lR_prod{k} * CB';
                R_prod{k}  = CB *  R_prod{k} * CB';

%               Rimuoviamo eventuale sporcizia numerica
                lR_prod{k} = pulisci_matrice(lR_prod{k}, 1e3 * eps);
                R_prod{k}  = pulisci_matrice( R_prod{k}, 1e3 * eps);

%               Scriviamo il superoperatore associato
                sop_comm_R{k} = genera_superoperatore_sx(lR_prod{k}) ...
                              - genera_superoperatore_dx(lR_prod{k}) ;
            end
        end
    end
end