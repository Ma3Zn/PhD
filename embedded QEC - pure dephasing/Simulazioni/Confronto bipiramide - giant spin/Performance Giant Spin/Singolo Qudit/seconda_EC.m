function errore_medio = seconda_EC(sim, rho, iter)
%   Funzione che dati in input i parametri della simulazione e la matrice
%   rho da stabilizzare applica un secondo ciclo incoerente di EC

%   Allocazione di p_osservazione ed errore. Gli indici si riferiscono al k
%   di partenza ed al k di proiezione
    p_osservazione = zeros(sim.dim_q/2,1);
    errore         = zeros(sim.dim_q/2,1);

%   Costruiamo la rho dell'ancilla
    psi_a    = zeros(sim.dim_a);
    psi_a(1) = 1;
    rho_a    = psi_a * psi_a';

%   Costruiamo la rho del sistema totale
    rho = kron(rho, rho_a);

%   Vettorizziamo la rho
    sim.rho = vettorizza_matrice(rho);

%   Simulazione gate logico (per poterlo fare devi calcolarti
%   l'opportuno stato finale, oppure simulare una sola fase di idle)
    sim.rho = sim.ev_gl * sim.rho;

%   Simulazione stabilizzazione
    sim.rho = sim.ev_CU * sim.rho;

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

%       Chiamata ricorsiva alla funzione per simulare la ripetizione di più
%       cicli di EC
        if (iter < sim.EC_iter)
            errore(k) = seconda_EC(sim, rho, iter + 1);
        else
%           Applichiamo una stabilizzazione ideale del qudit isolando di fatto
%           i soli errori NON correggibili introdotti dalla procedura di EC
            errore(k) = stabilizzazione_ideale(sim, rho);

% % % % %           Calcolo errore (correggibile + non correggibile) della procedura
% % % %             errore(k) = 1 - sim.psi_fin' * rho * sim.psi_fin;

% % % % %           Misura in base X logica
% % % %             rho = sim.HL * rho * sim.HL';
% % % %             rho = sim.RL * rho * sim.RL';
% % % %
% % % % %           Calcoliamo l'errore come l'errore che compiremmo eseguendo una 
% % % % %           misura logica del qudit al termine di un effettivo circuito
% % % %             errore(k) = calcola_errore_misura_logica_qudit(sim, rho);
        end

%       Assicuriamoci di non avere errori numerici che ci portino l'errore
%       ad essere negativo
        if (errore(k) < -1e1 * eps)
            errore(k) = abs(real(errore(k)));
            fprintf("\n\nerr negativo\n\n")
        end
    end

%   Calcoliamo l'errore medio
    errore_medio = p_osservazione.' * errore;

end