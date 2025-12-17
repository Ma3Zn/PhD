function errore_medio = ripeti_misura_ancilla(sim, rho, k_vot, iter)
%  Funzione ricorsiva per la ripetizione della misura dell'ancilla per
%  ottenere il k indicato per la recovery

    p_osservazione = zeros(1,sim.dim_a);
    errore = zeros(1,sim.dim_a);

%   Aggiorniamo l'indice iterativo di misura
    iter = iter + 1;

%   Tempo idle pre-misura
    rho_init = devettorizza_matrice(sim.ev_idle * vettorizza_matrice(rho));

    for k_gen = 1:sim.dim_a

%       Inizializzazione variabili
        k_vot_attuale = [k_vot, k_gen];

%       Probabilità di eseguire il k-esimo proiettore generalizzato
        p_osservazione(k_gen) = real(trace(sim.prj_gen{k_gen} * rho_init));

        if(p_osservazione(k_gen) == 0)
            continue;
        end

%       Inizializziamo le variabili
        p_prj = zeros(1,sim.dim_a);
        err = zeros(1,sim.dim_a);

        for k = 1:sim.dim_a

%           Inizializzazione variabili
            rho = rho_init;

%           Probabilità di eseguire il k-esimo proiettore e termine di
%           rinormalizzazione dello stato
            p_prj(k) = trace(sim.prj{k} * rho_init);

            if (~(p_prj(k) > 0))
                continue;
            end

%           Eseguiamo la proiezione del sistema
            rho = sim.prj{k} * rho * sim.prj{k}';

%           Normalizziamo la proiezione
            rho = rho / p_prj(k);

%           Verifichiamo se dobbiamo o meno eseguire un'altra misura
%           dell'ancilla
            if (iter < 4)
                err(k) = ripeti_misura_ancilla(sim, rho, k_vot_attuale, iter);
            else
%               Selezioniamo l'indice di recovery da applicare
                k_rec = seleziona_indice(k_vot_attuale);
    
%               Tracciamo fuori l'ancilla
                rho = partial_trace_1(rho, sim.dim_q, sim.dim_a);

%               Applichiamo la procedura di recovery incoerente al sistema
                rho = sim.ev_R{k_rec} * vettorizza_matrice(rho);
            
%               Devettorizziamo il risultato
                rho = devettorizza_matrice(rho);
    
%               Eseguiamo una tomografia logica del sistema
                err(k) = calcola_errore_tomografia_logica_qudit(sim, rho);
            end
        end

%       Calcoliamo l'errore per il k_gen attuale
        errore(k_gen) = p_prj * err.';
    end

%   Calcoliamo l'errore medio della procedura di ripetizione di misura
%   dell'ancilla, di fatto includendo dunque anche le probabilità di
%   osservazione dei vari k_gen
    errore_medio = p_osservazione * errore.';
end