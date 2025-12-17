% close all
clear all
clc

% Flag per il plot di overlap o errore -- [0, overlap] [altrimenti, errore]
% TODO:: far si di poter plottare entrambi su due figure diverse senza
% dover metterci il doppio del tempo necessario
flag_plot = 1;

% Definiamo la sequenza di spin che compongono il nostro anello/catena
N_spins  = [4 6 8 10 12];
N_spins  =  12;
periodic =  true;

% Numero di step di trotter da eseguire
N_it = 0:1:7;
% N_it = 0:5:35;
N_dt = length(N_it);

% Creiamo la figura su cui plottare i risultati
figure();
hold on
grid on

% Coupling Heisenberg
J.x = 1;
J.y = 1;
J.z = 1;

% velocità annealing
lambda = 0.1;

drive.f = @(t,par) (sin(par.lambda*t*pi/2))^(0.9);
% drive.f = @(t,par) par.lambda*t;

% Drive per il termine di splitting
drive_sp.f = @drive_splitting;

% Numero step overshooting
N_ov = 0;

% Esponente massimo della potenza dell'Hamiltoniana utilizzate per
% costruire lo spazio di krylov impiegato per la subpspace expansion
N_potenze = 2;

% Fattore di dilatazione temporale
for N_spin = N_spins

%   Teniamo traccia dell'avanzamento della simulazione
    N_spin

%   Generiamo parametri simulazione
    sim = genera_input(N_spin, periodic, J, drive, drive_sp);

%   Settiamo il parametro di velocità dell'annealing
    sim.drive.par.lambda = lambda;

%   Settiamo parametri drive ZFS
    sim.drive_sp.par.lambda = lambda;
    sim.drive_sp.par.n = N_spin;

%   Calcoliamo il ground state di H0 ed inizializziamo lì il sitema
    [V, e] = eig(sim.H0);
    
    psi_in = V(:,1);

%   Costo preparazione stato iniziale
    p_prep = 0;

%   Plottiamo gli stati iniziali
    dec2bin(find(V(:,1))-1)

%   Calcoliamo l'overlap al variare del numero di trotter step della
%   digitalizzazione del processo di annealing
    for i = 1:N_dt
    
        N = N_it(i);
    
%       Calcoliamo la profondità del circuito in funzione del numero di
%       step di trotterizzazione che dobbiamo eseguire. Sono 6 CNOT per
%       implementare il termine di Heisenberg di un trotter step.
        profondita(i) = (N + N_ov) * sim.trt.overhead + p_prep;
    
%       Calcoliamo l'evoluzione temporale del sistema usando la formula
%       di trotter del second'ordine
        tic
        psi_fin = calcola_evoluzione_trotter(N, sim, psi_in);
    
%       Eseguiamo overshooting
        psi_fin = overshoot(N, N_ov, sim, psi_fin);

%       Applichiamo la subspace expansion agli stati generati fin'ora
        psi_fin = subspace_expansion(sim, psi_fin, N_potenze);

%       Calcola overlap evoluzione
        [overlap(i), err(i)] = calcola_overlap(sim, psi_fin); 
        toc
    end

%   Verifichiamo quale plot eseguire
    if (flag_plot == 0)
%       Plottiamo l'overlap calcolato
        plot(profondita, overlap, 'LineWidth', 2.5);
    else
%       Plottiamo l'errore calcolato
        semilogy(profondita, err, 'LineWidth', 2.5);
    end
end

%%  plottiamo l'overlap in funzione della profondità del circuito

% Addobbiamo il plot
legend(string(N_spins), 'Location', 'southeast', 'FontSize', 16);

% Verifichiamo quale plot eseguire
if (flag_plot == 0)
    title(strcat('Overlap $|gs\rangle$ - ($\lambda=',string(sim.drive.par.lambda),'$) - $N = ', num2str(N_potenze), '$'), 'Interpreter', 'latex', 'FontSize',20);
    xlabel('Circuit Depth','Interpreter','latex',"FontSize",20); 
    ylabel('$|| \Pi_0|\psi\rangle ||_2$','Interpreter','latex',"FontSize",20,'Rotation',90);
    set(gca,'FontSize',15);
else
    title(strcat('Errore energia $|gs\rangle$ - ($\lambda=',string(sim.drive.par.lambda),'$) - $N = ', num2str(N_potenze), '$'), 'Interpreter', 'latex', 'FontSize',20);
    xlabel('Circuit Depth','Interpreter','latex',"FontSize",20);
    ylabel('$|| e_{gs} - e_{app} ||_2/||e_{gs}||_2$','Interpreter','latex',"FontSize",20,'Rotation',90);
    set(gca,'YScale','log','FontSize',15);
end