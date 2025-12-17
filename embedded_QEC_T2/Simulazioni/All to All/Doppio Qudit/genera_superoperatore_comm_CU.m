function sop_comm_CU = genera_superoperatore_comm_CU(dim, CB, q)
%   Funzione che dati in input la dimensione dei qudit e la matrice per il
%   cambiamento di base da Logica a computazionale, calcola il
%   superoperatore commutatore delle operazioni per la error detection

    I  = speye(dim);

%   costruiamo la forma base del CU
    CU = genera_CU(dim);
    switch q
        case 1
            CU_tot = kron(kron(CU, I), I);
        case 2
            CU_tot = kron(kron(I, CU), I);
        case 3
            CU_tot = kron(kron(I, I), CU);
        otherwise
            disp('ERRORE:: indice del qudit da correggere sbagliato');
    end

%   Calcoliamo il logaritmo di questo operatore
    lCU_tot = logm(full(CU_tot));

%   scriviamo il logaritmo in base computazionale
    lCU_tot = CB * sparse(lCU_tot) * CB';

%   Rimuoviamo eventuale sporcizia numerica
    tmp = pulisci_matrice(lCU_tot, 1e3 * eps);

%   Calcoliamo il superoperatore commutatore
    sop_comm_CU = genera_superoperatore_sx(tmp) ...
                - genera_superoperatore_dx(tmp) ;
end