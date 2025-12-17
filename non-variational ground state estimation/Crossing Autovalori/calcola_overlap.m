function [time, overlap] = calcola_overlap(sim, psi)
%   Funzione che calcola l'overlap del nostro stato finale con il ground
%   state effettivo dell'Hamiltoniana
    
%   ATTENZIONE:: L'overlap sull'intero sottospazio non va bene, ma fintanto
%   che hai un hamiltoniana non degenere questa funzione va bene ...

%   Parametri solver ODE
    Relative_error = 1e-8;
    Absolute_error = 1e-8;

%   Tempo finale simulazione
    tf = 1/sim.drive.par.lambda;
    delta_t = 1e-2 * tf;
    time = 0:delta_t:tf;
    Ntime = length(time);

    overlap = zeros(1,Ntime+1);

%   solution of the differential equations
    options = odeset('RelTol',Relative_error,'AbsTol',Absolute_error);
    [~, y] = ode45(@(T,y) evol(T, y, sim), time, psi, options);
    

%%  Calcolo overlap
    
    for it = 1:Ntime
        Ht = sim.H(sim.drive, time(it), sim.Hz, sim.Hxy) + sim.Hs;
        [Vt,Et] = eig(full(Ht));

 %      In caso ci sia un ground state degenere
        gs_e = min(diag(Et));
        idxs = find(diag(Et) - gs_e < 1e-12);

        for idx = idxs'
            overlap(it+1) = overlap(it+1) + abs(Vt(:,idx)' * y(it,:).')^2;
        end
    end
end