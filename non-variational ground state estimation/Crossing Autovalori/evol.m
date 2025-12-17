function dy = evol(T, y, sim)

%   Creazione dell'Hamiltoniana per l'evoluzione adiabatica del sistema
    Htot = sim.H(sim.drive, T, sim.Hz, sim.Hxy) + sim.Hs;

    dy = -1i * Htot * y;
end