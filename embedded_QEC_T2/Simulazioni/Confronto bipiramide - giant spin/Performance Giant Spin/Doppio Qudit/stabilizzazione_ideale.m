function errore = stabilizzazione_ideale(sim, rho_init)
%   Funzione che dati in input i parametri della simulazione e la rho
%   attuale esegue una stabilizzazione ideale dello stato del sistema 
%   ridotto e ritorna l'errore medio che osserviamo facendo ciò

%   Allocazione di p_osservazione ed errore. Gli indici si riferiscono al k
%   di partenza ed al k di proiezione
    err = zeros((sim.dim_q/2)^3, 1);
    p   = zeros((sim.dim_q/2)^3, 1);

    for kc = 1:sim.dim_a
        for kt = 1:sim.dim_a
            for ks = 1:sim.dim_a
%               Calcoliamo l'indice attuale
                k = (kc - 1)* sim.dim_a^2 + (kt-1) * sim.dim_a + ks;

%               Inizializziamo rho
                rho = rho_init;

%               Calcola probabilità osservazione
                p(k) = trace(sim.prj_rid{k} * rho);

%               Verifichiamo di avere una probbailità non nulla
                if ( ~(p(k) > 0) )
                    continue;
                end

%               Eseguiamo la proiezione del sistema
                rho = sim.prj_rid{k} * rho * sim.prj_rid{k}';
                
%               Normalizziamo lo stato del sistema
                rho = rho / p(k);

%               Applichiamo la procedura algebrica di recovery per avere un
%               semplice modo per il calcolo della fidelity finale
                rho = sim.R{k} * rho * sim.R{k}';

%               Calcolo errore
                err(k) = 1 - sim.psi_fin' * rho * sim.psi_fin;

%               Assicuriamoci di non avere errori numerici che ci portino
%               l'errore ad essere negativo
                if (err(k) < 0)
                    err(k) = abs(real(err(k)));
%                     fprintf("\n\nerr negativo\n\n");
                end
            end
        end
    end

    %   Calcoliamo e torniamo l'errore medio commesso
    errore = p.' * err;
end