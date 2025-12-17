% close all
clear all
clc

% Definiamo la sequenza di spin che compongono il nostro anello/catena
N_spins  = [8, 10];
N_spins  = 8;
periodic = true;

% Numero autovalori da plottare
N_aut = 5;

% Valori di lambda da testare
lambda_vals = 0.1;

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

%   Calcoliamo l'overlap al variare di lambda
    for lambda = lambda_vals
    
        sim.drive.par.lambda = lambda;
        sim.drive.par.n = N_spin;

        %   Tempo finale simulazione
        tf = 1/lambda;
        delta_t = 1e-2 * tf;
        time = 0:delta_t:tf;
        Ntime = length(time);

        autoval = zeros(N_aut, Ntime);

%       Includiamo il valore attuale di lambda all'interno dei parametri
%       della simulazione
        sim.drive.par.lambda = lambda;
    
%       Calcoliamo istante per istante per istante l'evoluzione degli
%       autovettori
        tic
        for it = 1:Ntime
            it
            tic
%           Calcoliamo l'Hamiltoniana attuale
            H = sim.H(sim, time(it));
            
%           Calcoliamo gli autovalori di H
            [~, e] = eig(full(H));

%           Riordiniamo gli autovalori 
            e = diag(e);

%           Recuperiamo gli N_aut autovalori più piccoli di H
            autoval(:,it) = e(1:N_aut);
            toc
        end
        toc

%%      PLOT AUTOVALORI

%       Creiamo la figura su cui plottare i risultati a lambda fissato
        figure()
        hold on
        grid on

%       Plottiamo l'evoluzione degli autovalori
        plot(time*lambda, autoval, 'LineWidth', 2.5);

%       legend(string(1:N_aut));
        xlabel('$t * \lambda$','Interpreter','latex',"FontSize",20);
        ylabel('$E$','Interpreter','latex',"FontSize",20,'Rotation',0);
        title(strcat('Eigenvalues($\lambda$) - $\lambda = ', string(lambda),'$ - $N_s =  ',string(N_spin), '$'), 'Interpreter', 'latex', 'FontSize',20);
    
% % % %%      PLOT DIFFERENZA PRIMI AUTOVALORI
% % % 
% % % %       Creiamo la figura su cui plottare i risultati a lambda fissato
% % %         figure()
% % %         hold on
% % %         grid on
% % % 
% % % %       Plottiamo l'evoluzione degli autovalori
% % %         plot(time*lambda, abs(autoval(1,:) - autoval(2,:)), 'LineWidth', 2.5);
% % % 
% % % %       legend(string(1:N_aut));
% % %         xlabel('$t * \lambda$','Interpreter','latex',"FontSize",20);
% % %         ylabel('$|E_0 - E_1|$','Interpreter','latex',"FontSize",20,'Rotation',90);
% % %         title(strcat('Eigenvalues($\lambda$) - $\lambda = ', string(lambda),'$ - $N_s =  ',string(N_spin), '$'), 'Interpreter', 'latex', 'FontSize',20);
    

    end
end