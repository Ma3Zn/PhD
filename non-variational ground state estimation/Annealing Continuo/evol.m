function dy = evol(T, y, sim)

%   Creazione dell'Hamiltoniana per l'evoluzione adiabatica del sistema
    Htot = sim.H(sim, T);

    dy = -1i * Htot * y;
end