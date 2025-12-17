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

    for k = 1:sim.dim_q/2

%       Devettorizziamo il nostro stato
        rho = devettorizza_matrice(sim.rho);

%       Calcola probabilità osservazione
        p_osservazione(k) = real(trace(sim.prj{k} * rho));

        if(p_osservazione(k) == 0)
            continue;
        end

%       Eseguiamo la proiezione del sistema
        rho = sim.prj{k} * rho * sim.prj{k}';

%       Normalizziamo la proiezione
        rho = rho / p_osservazione(k);

%       Tracciamo fuori l'ancilla
        rho = partial_trace_1(rho, sim.dim_q, sim.dim_a);

%       Applichiamo la procedura di recovery incoerente al sistema
        rho = sim.ev_R{k} * vettorizza_matrice(rho);
        
%       Devettorizziamo il risultato
        rho = devettorizza_matrice(rho);

% % % %       Eseguiamo un secondo ciclo di EC
% % %         errore(k) = seconda_EC(sim, rho, 1);

% % % %       Calcolo errore (correggibile + non correggibile) della procedura
% % %         errore(k) = 1 - sim.psi_fin' * rho * sim.psi_fin;

% % % %       Applichiamo una stabilizzazione ideale del qudit isolando di fatto
% % % %       i soli errori NON correggibili introdotti dalla procedura di EC
% % %         errore(k) = stabilizzazione_ideale(sim, rho);

%       Eseguiamo una tomografia logica del sistema
        errore(k) = calcola_errore_tomografia_logica_qudit(sim, rho);

% % % %       Misura in base X logica
% % %         rho = sim.HL * rho * sim.HL';
% % %         rho = sim.RL * rho * sim.RL';

% % % %       Calcoliamo l'errore come l'errore che compiremmo eseguendo una 
% % % %       misura logica del qudit al termine di un effettivo circuito
% % %         errore(k) = calcola_errore_misura_logica_qudit(sim, rho);

%       Assicuriamoci di non avere errori numerici che ci portino l'errore
%       ad essere negativo
        if (errore(k) < 0)
            errore(k) = abs(real(errore(k)));
%             fprintf("\n\nerr negativo\n\n")
        end

%       Tempo simulazioni
        T = T + sim.T_tot;
    end

    % tempo medio esecuzione simulazione
    T = T / (sim.dim/2);
end
