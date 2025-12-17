function sim = genera_input(N, periodic, J, drive, drive_sp)
%   Funzoine che genera tutti i parametri di input necessari per le
%   simulazioni che andremo ad eseguire

%   Costruiamo la catena di spin
    spin_chain = zeros(1,N) + 1/2;

%   Costruiamo le Hamiltoniane del nostro sistema
    [H0, Hs, H_heis] = costruisci_hamiltoniane(spin_chain, periodic);

%   Costruiamo la function handle per la valutazione dell'Hamiltoniana
    H = @(drive,t,J,Hx,Hy,Hz) J.z * Hz + drive.f(t,drive.par) * (J.x * Hx + J.y * Hy);

%%  Decomposizione Trotter

% % % %   First order trotter decomposition
% % %     trt.a = [1];
% % %     trt.b = [0];
% % %     trt.overhead = 3; % Numero CNOT per implementazione
% % %     trt.name = "Trotter (I)";

% % % %   Verlet o Leapfrog (II order -- Eff2 = 10.7) [meh]
% % %     trt.a = [1/2 1/2];
% % %     trt.b = [1 0];
% % %     trt.overhead = 9; % Numero CNOT per implementazione
% % %     trt.name = "Verlet (II)";

%   Omelyan (II order -- Eff2 = 29.2) [Ottima]
    trt.a = [0.1931833275037836 0.613633344992433 0.1931833275037836];
    trt.b = [1/2 1/2 0];
    trt.overhead = 15; % Numero CNOT per implementazione
    trt.name = "Omelyan (II)";

% % % %   Forest-Ruth (IV order -- Eff4 = 0.315) [Pessima]
% % %     trt.a = [0.6756035959798288 -0.175603595979829 -0.175603595979829 0.6756035959798288];
% % %     trt.b = [1.351207191959658 -1.702414383919316 1.351207191959658 0];
% % %     trt.overhead = 21; % Numero CNOT per implementazione
% % %     trt.name = "Forest-Ruth (IV)";

% % % %   Omelyan's Forest-Ruth-Type (IV order -- Eff4 = 4.24) [Ok]
% % %     trt.a = [0.1720865590295143 -0.1616217622107222 0.979070406362416 -0.1616217622107222 0.1720865590295143];
% % %     trt.b = [0.5915620307551568 -0.091562030755157 -0.091562030755157 0.5915620307551568 0];
% % %     trt.overhead = 27; % Numero CNOT per implementazione
% % %     trt.name = "Omelyan - FR (IV)";

% % % %   Suzuki (IV order -- Eff4 = 1.10) [Pessima]
% % %     trt.a = [0.2072453858971879 0.4144907717943757 -0.121736157691564 -0.121736157691564 0.4144907717943757 0.2072453858971879];
% % %     trt.b = [0.4144907717943757 0.4144907717943757 0.171018456411249 0.4144907717943757 0.4144907717943757 0];
% % %     trt.overhead = 33; % Numero CNOT per implementazione
% % %     trt.name = "Suzuki (IV)";

% % % %   Optimised 4th order (IV order -- Eff4 = 10.5) [Ottima lievemente out of budget]
% % %     trt.a = [0.09257547473195787 0.4627160310210738 -0.055291505753032 -0.055291505753032 0.4627160310210738 0.09257547473195787];
% % %     trt.b = [0.2540996315529392 -0.1676517240119692 0.827104184918060 -0.1676517240119692 0.2540996315529392 0];
% % %     trt.overhead = 33; % Numero CNOT per implementazione
% % %     trt.name = "Optimised (IV)";

% % % %   Blanes & Moan (IV order -- Eff4 = 10.2) [Ottima ma costosa]
% % %     trt.a = [0.07920369643119569 0.353172906049774 -0.0420650803577195 0.219376955753500 -0.0420650803577195 0.353172906049774 0.07920369643119569];
% % %     trt.b = [0.209515106613362 -0.143851773179818 0.434336666566456 0.434336666566456 -0.143851773179818 0.209515106613362 0];
% % %     trt.overhead = 39; % Numero CNOT per implementazione
% % %     trt.name = "BM (IV)";

% % % %   Blanes & Moan (VI order) [Ottima ma costosa]
% % %     trt.a = [0.0502627644003922 0.413514300428344 0.0450798897943977 -0.188054853819569 0.54196067845078 -0.725525558508690 0.54196067845078 -0.188054853819569 0.0450798897943977 0.413514300428344 0.0502627644003922];
% % %     trt.b = [0.148816447901042 -0.132385865767784 0.067307604692185 0.432666402578175 -0.016404589403618 -0.016404589403618 0.432666402578175 0.067307604692185 -0.132385865767784 0.148816447901042 0];
% % %     trt.overhead = 63; % Numero CNOT per implementazione
% % %     trt.name = "BM (VI)";

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
end