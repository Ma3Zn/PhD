function [CU, sop_comm_CU] = genera_superoperatore_comm_CU(dim, CB)
%   Funzione che dati in input la dimensione dei qudit e la matrice per il
%   cambiamento di base da Logica a computazionale, calcola il
%   superoperatore commutatore delle operazioni per la error detection

%   costruiamo la forma base del CU
    CU = genera_CU(dim);

%   Calcoliamo il logaritmo di questo operatore
    lCU = logm(CU);

%   scriviamo il logaritmo in base computazionale
    lCU = CB * lCU * CB';

%   Rimuoviamo eventuale sporcizia numerica
    lCU = pulisci_matrice(lCU, 1e3 * eps);

%   Calcoliamo il superoperatore commutatore
    sop_comm_CU = genera_superoperatore_sx(lCU) ...
                - genera_superoperatore_dx(lCU) ;
end