% close all
clear all
clc

% Definiamo la sequenza di spin che compongono il nostro anello/catena
N_spins  = [8, 10, 12];
N_spins  = 8;

periodic = true;

% Valori di lambda da testare
lambda_vals = logspace(0,-2,5);

% Coupling Heisenberg
J = 1;

drive.sp = @(t,par) ( tan(1) - tan((par.lambda*t)^par.n) ); 
drive.f  = @(t,par) ( sin(par.lambda*t*pi/2))^(0.9);

% Ciclo sulle diverse catene di spin
for N_spin = N_spins

    N_spin

%   Controlliamo che N_spin sia pari
    if (mod(N_spin,2)==1)
        continue;
    end

%   Generiamo parametri simulazione
    sim = genera_input(N_spin, periodic, J, drive);

    sim.drive.par.lambda = 0;
    sim.drive.par.n      = N_spin;

%   Generiamo lo stato iniziale del sistema
    sim.psi_in = genera_stato_iniziale(sim);

%   Creiamo la figura su cui plottare i risultati al variare di lambda
    figure()
    hold on
    grid on

%   Calcoliamo l'overlap al variare di lambda
    for lambda = lambda_vals
    
        sim.drive.par.lambda = lambda;

%       Calcola overlap evoluzione
        tic
        [time, overlap] = calcola_overlap(sim);
        toc

%       plottiamo l'overlap attuale
        plot(time * lambda, overlap, 'LineWidth', 2.5);

    end

%   Addobbiamo il plot
    legend(string(lambda_vals), 'Location', 'southwest', 'FontSize', 16);
    title(strcat('Overlap($\lambda$) - Ground State'), 'Interpreter', 'latex', 'FontSize',20);
    xlabel('$t * \lambda$','Interpreter','latex',"FontSize",20);
    ylabel('$|<\psi|\mathrm{gs}>|$','Interpreter','latex',"FontSize",20,'Rotation',90);
end