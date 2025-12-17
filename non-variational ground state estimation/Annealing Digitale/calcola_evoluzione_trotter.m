function psi = calcola_evoluzione_trotter(N, sim, psi)
%   Calcoliamo l'evoluzione temporale del sistema usando la formula di
%   trotter del second'ordine. N è il numero di step da dover eseguire.

    if (N == 0)
        return;
    end

%   Calcoliamo dt
    tf =  1 / sim.drive.par.lambda;
    dt = tf / N;
    time = 0:dt:tf-dt;
    Ntime = length(time);

%   Variabili d'appoggio per migliorare la leggibilità
    H_heis = sim.H_heis;
    J = sim.J;

%   Eseguiamo l'evoluzione del sistema
    for it = 1:Ntime

%       Tempo simulazione
        t = time(it);

%       Calcoliamo il centro di valutazione temporale
        t_eval = t+dt/2;

%       Calcoliamo l'Hamiltoniana
        H_odd  = sim.H(sim.drive, t_eval, J, H_heis.Sx_odd , H_heis.Sy_odd , H_heis.Sz_odd );
        H_even = sim.H(sim.drive, t_eval, J, H_heis.Sx_even, H_heis.Sy_even, H_heis.Sz_even);

%       Costruiamo gli operatori di evoluzione temporale come nella formula
%       di trotter al second'ordine
        if (mod(it, 2) == 0)
%           Calcoliamo l'operatore di evoluzione dato dalla
%           trotterizzazione opportuna
            U = calcola_decomposizione_trotter(sim.trt, H_odd, H_even, dt);
        else
%           Calcoliamo l'operatore di evoluzione dato dalla
%           trotterizzazione opportuna
            U = calcola_decomposizione_trotter(sim.trt, H_even, H_odd, dt);
        end

%       Calcoliamo l'operatore a singolo corpo di splitting per l'hamiltoninaa
%       iniziale che rimuoviamo da quella target
        Uz = fastExpm(-1i * sim.drive_sp.f(t_eval, sim.drive_sp.par) * sim.Hs * dt/2);

%       Aggiungiamo le rotazioni dovute alla parte a singolo corpo
        U = Uz * U * Uz;

%       Evolviamo il sistema con la formula di Trotter al second'ordine
        psi = U * psi; 
    end
end