function [T_rabi, T_esecuzione] = esegui_protocollo_sensing(sim, lv, T_max)
%   Funzione che esegue il protocollo di sensing sul sistema passato in
%   input ne ritorna il campo misurato

%   indici rho vettorizzata delle popolazioni di |0,k> per ogni k
    supp0 = sim.supp0_v;

%   indici rho vettorizzata delle popolazioni di |0,k> per ogni k
    supp1 = sim.supp1_v;

%   Prima rotazione
    sim.rho = sim.sop_rot * sim.rho;

%   Aggiorniamo il contatore che tiene traccia del numero di procedure di
%   QEC eseguite
    sim.it = sim.it + 1;

%   Eseguiamo un ciclo di QEC
    sim.rho = esegui_QEC(sim);

%   Ciclo principale di evoluzione del sistema in cui ci avvicianiamo al
%   minimo tramite ripetizioni di rot(dt) + EC
    while( real(sum(sim.rho(supp0))) > lv )
  
%       Controlliamo di non star divergendo, i.e., se NON c'è crossing
        if (sim.T > T_max || sim.T == inf)
            T_esecuzione = 0;
            T_rabi = 0;
            return;
        end

% % % %       Teniamo traccia della durata della rotazione logica CON EC
% % %         sim.T = sim.T + sim.dt + 2 * sim.durata_op_EC + sim.durata_misura;

%       Teniamo traccia della durata della rotazione logica SENZA EC
        sim.T = sim.T + sim.dt;

%       Aggiorniamo la rho_old
        sim.rho_old = sim.rho;

%       Eseguiamo un'ulteriore pezzo di rotazione
        sim.rho = sim.sop_rot * sim.rho;

%       Aggiorniamo il contatore che tiene traccia del numero di procedure
%       di QEC eseguite
        sim.it = sim.it + 1;

%       Eseguiamo un ciclo di QEC
        sim.rho = esegui_QEC(sim);

    end

%   Calcoliamo il tempo totale dell'esperimento (Tempo di rotazione logica
%   + tempo durata QEC + tempo durata misura QEC)
% % %     T_esecuzione = sim.T;
    T_esecuzione = sim.T + sim.it * (2 * sim.durata_op_EC + sim.durata_misura);

% % % %   __TPAR__
% % % %   Controlliamo quale definizione di durata utilizzare
% % %     if(~(sim.dt > 0))
% % %         sim.T = sim.it * sim.durata_op_EC;
% % %     end

% % % %   __TTOT__
% % % %   Aggiungiamo la durata degli step di EC al tempo di evoluzione
% % %     sim.T = sim.T + sim.it * sim.durata_op_EC;
 
%   Siamo usciti dal ciclo e dunque invertiamo il coseno per ricavare la
%   frequenza effettiva da cui poi calcolare il valore di Bx, il 2*pi ci va
%   perchè acos ritorna un valore in radianti
    T_rabi = sim.T * 2 * pi / acos( (real(sum(sim.rho(supp0))) - 1/2) * 2);
end