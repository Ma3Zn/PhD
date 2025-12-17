% close all
clear
clc

% Definiamo la sequenza di spin che compongono il nostro anello/catena
N_spin  =  8;
periodic =  true;

% Numero di step di trotter da eseguire
N_it = 3:1:7;
N_dt = length(N_it);

% Coupling Heisenberg
J.x = 1;
J.y = 1;
J.z = 1;

% velocità annealing
lambda = 0.1;

% drive.f = @(t,par) (sin(par.lambda*t*pi/2))^(0.9);
drive.f = @(t,par) par.lambda*t;

% Drive per il termine di splitting
drive_sp.f = @drive_splitting;

% Esponente massimo della potenza dell'Hamiltoniana utilizzate per
% costruire lo spazio di krylov impiegato per la subpspace expansion
N_potenze = 1;

% Selezioniamo quali valori di aspettazione caricare
dir = "Data/";

data = "25-03-27";

backend_noise = "/Fez";
% backend_noise = "/Future";

dir = strcat(dir, data, backend_noise);

% Ordine della SbE da utilizzare
ordine_SbE = 1;

% Indice di filtro per i folding factors
data_range = 1:3:17;

% Precisione stima valori di aspettazione misurati
prec_est = 1e-2;

% Numero di sample da generare per ogni campionamento dei valori di
% aspettazione noisy
N_campioni_noisy = 1e3;

% Generiamo parametri simulazione
sim = genera_input(N_spin, periodic, J, drive, drive_sp);

% Settiamo il parametro di velocità dell'annealing
sim.drive.par.lambda = lambda;

% Settiamo parametri drive ZFS
sim.drive_sp.par.lambda = lambda;
sim.drive_sp.par.n = N_spin;

% Calcoliamo l'energia del ground state effettivo
sim.e_gs = calcola_energia_gs(sim);

% Calcoliamo l'overlap al variare del numero di trotter step della
% digitalizzazione del processo di annealing
for i = 1:N_dt

%   Teniamo traccia dell'avanzamento della simulazione
    N_tstep = N_it(i)

%   Calcoliamo la profondità del circuito in funzione del numero di
%   step di trotterizzazione che dobbiamo eseguire.
    profondita(i) = N_tstep * sim.trt.overhead;

%   Generiamo l'imput per la SbE relativa all'attuale numero di Trotter
%   Step
    SbE = genera_input_SbE(dir, N_tstep, ordine_SbE, data_range, prec_est, ...
                           N_campioni_noisy);

%   Calcoliamo le estrapolazioni dei vari valori di aspettazione
    exp_val_ZNE = calcola_valori_aspettazione_ZNE(SbE.ZNE.exp_val, SbE.ZNE.ff);

%   Costruiamo le matrici della SbE tramite ZNE
    [HH, SS] = costruisci_matrici_SbE_ZNE(exp_val_ZNE);

%   Calcoliamo la stima dell'energia tramite SbE ZNE
    e_SbE = esegui_SbE(HH, SS);

%   Calcola overlap evoluzione
%   TODO:: Spostare fuori dal loop la diagonalizzazione dell'Hamiltoniana
    err(i) = norm(e_SbE - sim.e_gs) / norm(sim.e_gs);
end

%% Plot

% Creiamo la figura su cui plottare i risultati
figure();
hold on
grid on

% Plottiamo l'errore calcolato
semilogy(profondita, err, 'LineWidth', 2.5);

%%  plottiamo l'overlap in funzione della profondità del circuito

% Addobbiamo il plot
legend(string(N_spin), 'Location', 'southeast', 'FontSize', 16);

title(strcat('Errore energia $|gs\rangle$ - ($\lambda=',string(sim.drive.par.lambda),'$) - $N = ', num2str(N_potenze), '$'), 'Interpreter', 'latex', 'FontSize',20);
xlabel('Circuit Depth','Interpreter','latex',"FontSize",20);
ylabel('$|| e_{gs} - e_{app} ||_2/||e_{gs}||_2$','Interpreter','latex',"FontSize",20,'Rotation',90);
set(gca,'YScale','log','FontSize',15);
