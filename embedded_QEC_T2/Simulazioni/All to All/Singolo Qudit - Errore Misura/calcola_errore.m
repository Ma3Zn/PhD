function [p_osservazione, errore, T] = calcola_errore(sim)
%   Funzione che forniti in input uno spin S ed un valore di T2
%   simula l'eovoluziione temporale dello spin S implementando una
%   decomposizione in rotazioni piane e Z delle varie unitarie che
%   compongono il circuito
%
%   Fornisce in output le performance del codice, i.e., le probabilità di
%   proiettare in un dato k e l'errore commesso per tale proiezione

%   Allocazione di p_osservazione ed errore. Gli indici si riferiscono al k
%   di partenza ed al k di proiezione
    p_osservazione = zeros(sim.dim_q/2,1);
    errore         = zeros(sim.dim_q/2,1);

%   Tempo simulaizoni
    T = 0;

%   Simulazione gate logico
    sim.rho = sim.ev_gl * sim.rho;

%   Simulazione stabilizzazione
    sim.rho = sim.ev_CU * sim.rho;

%   Tempo idle pre-misura
    sim.rho = sim.ev_idle * sim.rho;

    rho_init = devettorizza_matrice(sim.rho);

    for k_gen = 1:sim.dim_a

%       Probabilità di eseguire il k-esimo proiettore generalizzato
        p_osservazione(k_gen) = real(trace(sim.prj_gen{k_gen} * rho_init));

        if(p_osservazione(k_gen) == 0)
            continue;
        end

%       Inizializziamo le variabili
        p_prj = zeros(1,sim.dim_a);
        err = zeros(1,sim.dim_a);

        for k = 1:sim.dim_a
%           Probabilità di eseguire il k-esimo proiettore e termine di
%           rinormalizzazione dello stato
            p_prj(k) = trace(sim.prj{k} * rho_init);

            if (~(p_prj(k) > 0))
                continue;
            end

%           Eseguiamo la proiezione del sistema
            rho = sim.prj{k} * rho_init * sim.prj{k}';

%           Normalizziamo la proiezione
            rho = rho / p_prj(k);

            % DBG
%             trace(rho) - 1

%           ATTENZIONE!!

%           Ripetiamo la misura dell'ancilla per ottenere l'errore medio
%           che commettiamo con la proiezione in questo k
            err(k) = ripeti_misura_ancilla(sim, rho, k_gen, 1);
 
% % % %           Per NON ripetere la misura generalizzata dell'ancilla
% % % %           Selezioniamo l'indice di recovery da applicare
% % %             k_rec = k_gen;
% % % 
% % % %           Tracciamo fuori l'ancilla
% % %             rho = partial_trace_1(rho, sim.dim_q, sim.dim_a);
% % % 
% % % %           Applichiamo la procedura di recovery incoerente al sistema
% % %             rho = sim.ev_R{k_rec} * vettorizza_matrice(rho);
% % %         
% % % %           Devettorizziamo il risultato
% % %             rho = devettorizza_matrice(rho);
% % % 
% % % %           Eseguiamo una tomografia logica del sistema
% % %             err(k) = calcola_errore_tomografia_logica_qudit(sim, rho);
        end

%       Calcoliamo l'errore per il k_gen attuale
        errore(k_gen) = p_prj * err.';

%       Assicuriamoci di non avere errori numerici che ci portino l'errore
%       ad essere negativo
        if (errore(k_gen) < 0)
            errore(k_gen) = abs(real(errore(k)));
%           fprintf("\n\nerr negativo\n\n")
        end 
    end
end
