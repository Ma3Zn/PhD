function psi = overshoot(N_trotter, N_overshoot, sim, psi)
%   Funzione che esegue l'overshooting del processo di anneling

    if(N_trotter == 0)
        return;
    end

%   Calcoliamo il dt attuale
    t_fin = 1 / sim.drive.par.lambda;
    dt = t_fin / (N_trotter * 1e0);

%   Variabili d'appoggio per migliorare la leggibilità
    H_heis = sim.H_heis;
    J = sim.J;

%   Costruiamo le Hamiltoniane da trotterizzare all'istante finale
    H_odd  = J.z*H_heis.Sz_odd  + J.x*H_heis.Sx_odd  + J.y*H_heis.Sy_odd;
    H_even = J.z*H_heis.Sz_even + J.x*H_heis.Sx_even + J.y*H_heis.Sy_even;

%   Calcoliamo l'operatore tempo indipendente che rimuove la degenerazione
%   del ground state
    Uz = fastExpm(-1i * sim.Hs * dt/2);

%   Tempo iniziale overshooting
    t = t_fin;

%   Cicliamo sul numero di step di overshooting eseguire
    for n = 1:N_overshoot
%       Istante temporale in cui valutare l'overshooting
        t_eval = t + dt/2;

%       Valutiamo le Hamiltoninane
        H1 = sim.drive.f(t_eval, sim.drive.par) * H_odd ;
        H2 = sim.drive.f(t_eval, sim.drive.par) * H_even;

        if (mod(n, 2) == 0)
            U1 = fastExpm(-1i * H1 * dt/2);
            U2 = fastExpm(-1i * H2 * dt  );

%           Calcoliamo l'operatore di evoluzione totale
            U = Uz * U1 * U2 * U1 * Uz;
        else
            U1 = fastExpm(-1i * H1 * dt  );
            U2 = fastExpm(-1i * H2 * dt/2);

%           Calcoliamo l'operatore di evoluzione totale
            U = Uz * U2 * U1 * U2 * Uz;
        end
%       Evolviamo il sistema con la formula di Trotter al second'ordine
        psi = U * psi; 
    end
end