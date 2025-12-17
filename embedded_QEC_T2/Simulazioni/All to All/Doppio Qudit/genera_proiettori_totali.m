function prj_tot = genera_proiettori_totali(dim, CB, q)
%   Funzione che dati in input la dimesione degli oggetti logici e la
%   matrice del cambio di base [Logica --> Computazionele], genera i
%   proiettori per la misyra delle varie ancille

%   Generiamo gli operatori di proiezione per un unico sistema
%   qudit-ancilla
    prj_a = genera_proiettori_ancilla(dim);
    prj_i = genera_proiettori_ideali(dim);

%   Allochiamo lo spazio opportuno per i proiettori totali
    prj_tot = cell(1,(dim/2)^3);

    for kc = 1 : dim/2
        for kt = 1 : dim/2
            for ks = 1 : dim/2
                k = (kc -1) * dim^2/4 + (kt - 1) * dim/2 + ks;

%               Costruiamo il prodotto dei vari operatori in base
%               all'indice del qudit su cui effetturare la procedura
%               effettiva di EC
                switch q
                    case 1
                        prj_tot{k} = kron(prj_a{kc}, kron(prj_i{kt}, prj_i{ks}));
                    case 2
                        prj_tot{k} = kron(prj_i{kc}, kron(prj_a{kt}, prj_i{ks}));
                    case 3
                        prj_tot{k} = kron(prj_i{kc}, kron(prj_i{kt}, prj_a{ks}));
                    otherwise
                        disp('ERRORE:: indice del qudit da correggere sbagliato');
                end

%               riscriviamo gli operatori in base computazionale
                prj_tot{k}  = CB * prj_tot{k}  * CB';

%               Rimuoviamo eventuale sporcizia numerica
                prj_tot{k}  = pulisci_matrice(prj_tot{k}, 1e3 * eps);

            end
        end
    end
end