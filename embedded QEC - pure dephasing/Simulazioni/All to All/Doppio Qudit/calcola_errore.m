function [p_osservazione, errore] = calcola_errore(sim, idx_q)
%   Funzione che forniti in input uno spin S ed un valore di T2
%   simula l'eovoluziione temporale dello spin S implementando una
%   decomposizione in rotazioni piane e Z delle varie unitarie che
%   compongono il circuito
%
%   Fornisce in output le performance del codice, i.e., le probabilità di
%   proiettare in un dato k e l'errore commesso per tale proiezione

%   Allocazione di p_osservazione ed errore. Gli indici si riferiscono al k
%   di partenza ed al k di proiezione
    p_osservazione = zeros((sim.dim_q/2)^3, 1);
    errore         = zeros((sim.dim_q/2)^3, 1);

%   Simulazione gate logico
    sim.rho = sim.sop_evol_gl * sim.rho;
    sim.rho = sim.sop_evol_gl * sim.rho;

%   Simulazione CU
    sim.rho = sim.sop_evol_CU * sim.rho;

%   Simuliamo la misura del sistema [TBD]
%   Va modificata in quanto deve prevedere la proiezione delle ancille,
%   una loro seguente traccia, la recovery incoerente e lo step di EC
%   ideale del sistema
    for k = 1:sim.dim_a^3

%       Devettorizziamo il nostro stato
        rho = devettorizza_matrice(sim.rho);

%       Calcola probabilità osservazione
        p_osservazione(k) = trace(sim.prj_tot{k} * rho);

        if(~(p_osservazione(k) > 0))
            continue;
        end

%       Eseguiamo la proiezione del sistema
        rho = sim.prj_tot{k} * rho * sim.prj_tot{k}';

%       Normalizziamo la proiezione
        rho = rho / p_osservazione(k);

%       Tracciamo fuori le ancille
        rho_rid = traccia_ancille(sim, rho, idx_q);

%       Applichiamo la procedura di recovery incoerente al sistema
        rho_rid = sim.sop_evol_R{k} * vettorizza_matrice(rho_rid);

%       Devettorizziamo il risultato
        rho_rid = devettorizza_matrice(rho_rid);

%       Applichiamo una stabilizzazione ideale del qudit isolando
%       di fatto i soli errori NON correggibili introdotti dalla
%       procedura di EC
        errore(k) = stabilizzazione_ideale(sim, rho_rid);
    end
end