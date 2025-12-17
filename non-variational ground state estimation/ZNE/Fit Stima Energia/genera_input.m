function sim = genera_input(N, periodic, J, drive, drive_sp)
%   Funzoine che genera tutti i parametri di input necessari per le
%   simulazioni che andremo ad eseguire

%   Costruiamo la catena di spin
    spin_chain = zeros(1,N) + 1/2;

%   Costruiamo le Hamiltoniane del nostro sistema
    [H0, Hs, H_heis] = costruisci_hamiltoniane(spin_chain, periodic);

%   Costruiamo la function handle per la valutazione dell'Hamiltoniana
    H = @(drive,t,J,Hx,Hy,Hz) J.z * Hz + drive.f(t,drive.par) * (J.x * Hx + J.y * Hy);

%%  Formula di Trotter Utilizzata

%   Omelyan (II order -- Eff2 = 29.2) [Ottima]
    trt.a = [0.1931833275037836 0.613633344992433 0.1931833275037836];
    trt.b = [1/2 1/2 0];
    trt.overhead = 15; % Numero CNOT per implementazione
    trt.name = "Omelyan (II)";

%%  Inizializzazione Parametri

%   Riempiamo la struttura di output con tutti i dati necessari
    sim.n               = N;
    sim.H0              = H0;
    sim.Hs              = Hs;
    sim.H_heis          = H_heis;
    sim.H               = H;
    sim.drive           = drive;
    sim.drive_sp        = drive_sp;
    sim.J               = J;
    sim.trt             = trt;
end