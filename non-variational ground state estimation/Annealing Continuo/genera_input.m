function sim = genera_input(N, periodic, J, drive)
%   Funzoine che genera tutti i parametri di input necessari per le
%   simulazioni che andremo ad eseguire

%   Costruiamo le Hamiltoniane del nostro sistema
    [Hs, Hz, Hxy] = costruisci_hamiltoniane(N, J, periodic);

%   Costruiamo la function handle per la valutazione dell'Hamiltoniana
    H = @(sim,t) sim.Hz + sim.drive.f(t,sim.drive.par) * sim.Hxy ...
               + sim.drive.sp(t,sim.drive.par) * sim.Hs;

%   Riempiamo la struttura di output con tutti i dati necessari
    sim.n        = N;
    sim.Hs       = Hs;
    sim.Hz       = Hz;
    sim.Hxy      = Hxy;
    sim.H        = H;
    sim.drive    = drive;
end