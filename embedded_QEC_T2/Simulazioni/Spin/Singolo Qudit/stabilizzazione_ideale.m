function errore = stabilizzazione_ideale(sim, rho)
%   Funzione che dati in input i parametri della simulazione e la rho
%   attuale esegue una stabilizzazione ideale dello stato del sistema e
%   ritorna l'errore medio che osserviamo facendo ciò

%   Recuperiamo il numero di osservazioni da esequire
    N_k = sim.dim_q/2;

%   Allochiamo lo spazio per errori e rispettive probabilità
    err = zeros(N_k, 1);
    p   = zeros(N_k, 1);   

%   Salviamo la rho iniziale
    rho_init = rho;

%   Cicliamo sui vari valori di k possibili
    for k = 1:N_k
%       Inizializiamo rho
        rho = rho_init;

%       Calcoliamo la probbailità di osservazione
        p(k) = trace(rho * sim.prj_q{k});

%       Verifichiamo di avere una probbailità non nulla
        if ( ~(p(k) > 0) )
            continue;
        end
        
%       Proiettiamo il sistema
        rho = sim.prj_q{k} * rho * sim.prj_q{k};

%       Normalizziamo lo stato del sistema
        rho = rho / p(k);

%       Applichiamo l'opportuno operatore di recovery per avere un semplice
%       modo di calcolare l'errore
        rho =  sim.R{k} * rho * sim.R{k}';

%       Calcoliamo l'errore commesso
        err(k) = 1 - sim.psi_fin' * rho * sim.psi_fin;

%       Assicuriamoci di non avere errori numerici che ci portino l'errore
%       ad essere negativo
        if (err(k) < -1e1 * eps)
            err(k) = abs(real(err(k)));
%             fprintf("\n\nerr negativo\n\n")
        end
    end

%   Calcoliamo e torniamo l'errore medio commesso
    errore = p.' * err;
end