function prj_rid = genera_proiettori_ridotti(dim, CB)
%   Funzione che dati in input la dimesione degli oggetti logici e la
%   matrice del cambio di base [Logica --> Computazionele], genera i
%   proiettori per la stabilizzazione ideale dei vari qudit

%   Generiamo gli operatori di proiezione per un unico qudit
    prj = genera_proiettori_ideali(dim);

%   Allochiamo lo spazio opportuno per i proiettori ridotti
    prj_rid = cell(1,(dim/2)^3);

    for kc = 1 : dim/2
        for kt = 1 : dim/2
            for ks = 1 : dim/2
                k = (kc -1) * dim^2/4 + (kt - 1) * dim/2 + ks;

%               Costruiamo il prodotto dei vari operatori
                prj_rid{k} = kron(prj{kc}, kron(prj{kt}, prj{ks}));

%               riscriviamo gli operatori in base computazionale
                prj_rid{k}  = CB * prj_rid{k}  * CB';

%               Rimuoviamo eventuale sporcizia numerica
                prj_rid{k}  = pulisci_matrice(prj_rid{k}, 1e3 * eps);

            end
        end
    end
end