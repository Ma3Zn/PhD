%% Pulizia ambiente di lavoro
close all
clear
clc

warning('off','MATLAB:legend:IgnoringExtraEntries');
warning('off','MATLAB:plot:IgnoreImaginaryXYPart');
warning('off','MATLAB:logm:nonPosRealEig');

%% Settiamo i parametri del grafico

t_opt = [5e-2 5e-6 5e-10 5e-14];
% t_opt = [3e-2 5e-2];

T2_vals = [1e5, 1e4, 1e3, 1e2];
T2_vals = 1e8;

% Valori angoli per un set di rotazioni piane universali (singolo corpo)
angle_vals = [pi/4 pi; pi/2 pi; pi/2 -pi/2; pi/2 -pi/4; pi/2 -pi/8]';
angle_vals = [0 0]';

% Possibili stati iniziali del qubit logico
init_superposition = [1 0; 0 1; 1/sqrt(2) 1/sqrt(2); 1/sqrt(2) -1/sqrt(2); 1/sqrt(2) 1i/sqrt(2); 1/sqrt(2) -1i/sqrt(2)]';
init_superposition = [1/sqrt(2) -1i/sqrt(2)]';

% Range di spin di cui fare il grafico
dim_vals = 4:2:12;

%% per generare correttamente l'input del programma
genera_input_1_2;

%% Facciamo i grafici della procedura di EC per tutti gli spin da testare
errore_medio = zeros(length(dim_vals),1);

for T2 = T2_vals
    T2

    fig = figure();

    for t = t_opt
        idx1 = 1;

        tic
        for dim = dim_vals
            dim
            errore_medio(idx1) = grafico_errore_medio(dim, T2, t, angle_vals, init_superposition);
            idx1 = idx1 + 1;
        end
        toc
    

%         errore_medio = pulisci_matrice(errore_medio, 1e3*eps);
        semilogy(dim_vals, errore_medio, '-x', 'MarkerSize', 10, 'LineWidth', 3);
        grid on
        hold on
    end

    % Salviamo il grafico
    title(strcat('$T_2 = ',num2str(T2),'$ns'),'Interpreter','latex');
    legend('$t = 5\times10^{-2}$', '$t = 5\times10^{-6}$', '$t = 5\times10^{-10}$', '$t = 5\times10^{-14}$', 'Interpreter','latex','Location', 'north','FontSize',16.2);
    xticklabels({'4', '', '6', '', '8', '', '10', '', '12','','14','','16'});
    xlabel('\textit{lv}','Interpreter','latex','FontSize',24);
    ylabel('$\bar{\mathcal{E}}_e$','Interpreter','latex','FontSize',24,'Rotation',0);
    
    saveas(fig, strcat('./images/pt_ef_andamento_S_',num2str(T2),'.fig'));
end