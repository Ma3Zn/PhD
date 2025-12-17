function sim = genera_input(N, periodic, J, drive, drive_sp, data)
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

%%  Recuperiamo i valori di aspettazione per la ZNE

    load(strcat(data, "/exp_values_4_tstep.txt"));
    load(strcat(data, "/exp_values_5_tstep.txt"));
    load(strcat(data, "/exp_values_6_tstep.txt"));
    load(strcat(data, "/folding_factors.txt"));

    data_range = 1:18;
    data_range = [1:5, 10:2:18];

    exp_values_4_tstep = exp_values_4_tstep(:, data_range);
    exp_values_5_tstep = exp_values_5_tstep(:, data_range);
    exp_values_6_tstep = exp_values_6_tstep(:, data_range);
    folding_factors    = folding_factors(data_range);

    % Riordiniamo i dati in modo opportuno
    [ZNE.ff, p]    = sort(folding_factors);
    ZNE.exp_val{1} = exp_values_4_tstep(:, p);
    ZNE.exp_val{2} = exp_values_5_tstep(:, p);
    ZNE.exp_val{3} = exp_values_6_tstep(:, p);

%%  Inizializzazione Parametri

%   Riempiamo la struttura di output con tutti i dati necessari
    sim.n        = N;
    sim.H0       = H0;
    sim.Hs       = Hs;
    sim.H_heis   = H_heis;
    sim.H        = H;
    sim.drive    = drive;
    sim.drive_sp = drive_sp;
    sim.J        = J;
    sim.trt      = trt;
    sim.ZNE      = ZNE;
end