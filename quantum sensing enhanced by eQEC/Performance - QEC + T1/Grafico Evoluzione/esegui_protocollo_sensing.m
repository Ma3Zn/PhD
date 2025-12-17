function sim = esegui_protocollo_sensing(sim)
%   Funzione che esegue il protocollo di sensing sul sistema passato in
%   input ne ritorna il campo misurato

%   indici rho vettorizzata delle popolazioni di |0,k> per ogni k
    supp0 = sim.supp0_v;

%   Inizializziamo il vettore contenente le popolazioni di |0>
    sim.pop(1) = 1;

%   Contatore
    sim.it = 2;

%   Ciclo principale di evoluzione del sistema in cui ci avvicianiamo al
%   minimo tramite ripetizioni di rot(dt) + EC
    while( sim.T < sim.T_fin )

%       Eseguiamo un'ulteriore pezzo di rotazione
        sim.rho = sim.sop_rot * sim.rho;

%       Eseguiamo un ciclo di QEC
        sim.rho = esegui_QEC(sim);

%       Inseriamo la popolazione logica di |0> dopo la prima rotazione
        sim.pop(sim.it) = real(sum(sim.rho(supp0)));

%       Teniamo traccia della durata dell'evoluzione del sistema
        sim.T = sim.T + sim.dt_tot;

%       Aggiorniamo il contatore che tiene traccia del numero di procedure
%       di QEC eseguite
        sim.it = sim.it + 1;
    end
end