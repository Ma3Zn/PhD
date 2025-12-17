%% Pulizia ambiente di lavoro
clear
clc

warning('off','MATLAB:legend:IgnoringExtraEntries');
warning('off','MATLAB:plot:IgnoreImaginaryXYPart');
warning('off','MATLAB:logm:nonPosRealEig');

%% Settiamo i parametri del grafico

% Apriamo la figura su cui fare il grafico
fig = figure(1);

% T2_vals = [1e7, 1e6, 1e5, 1e4, 1e3, 1e2];
T2_vals = [1e5, 1e4, 1e3, 1e2];
% T2_vals = 1e4

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

    idx1 = 1;
    
    tic
    for dim = dim_vals
        dim
        errore_medio(idx1) = grafico_errore_medio(dim, T2, angle_vals, init_superposition);
        idx1 = idx1 + 1;
    end
    toc
    
%     errore_medio = pulisci_matrice(errore_medio, 1e3*eps);
    semilogy(dim_vals, errore_medio, '-x', 'MarkerSize', 10, 'LineWidth', 3);
    grid on
    hold on
end

%%

% Salviamo il grafico
legend('$T_2 = 100\mu s$', '$T_2 = 10\mu s$', '$T_2 = 1 \mu s$', '$T_2 = 0.1\mu s$', 'Interpreter','latex','Location', 'southwest','FontSize',16.2);
% legend('T2 = 1e7', 'T2 = 1e6', 'T2 = 1e5', 'T2 = 1e4', 'T2 = 1e3', 'T2 = 1e2', 'Location', 'southwest');
xticklabels({'4', '', '6', '', '8', '', '10', '', '12'});
xlabel('\textit{lv}','Interpreter','latex','FontSize',24);
ylabel('$\bar{\mathcal{E}}_e$','Interpreter','latex','FontSize',24,'Rotation',0);

% saveas(fig, './images/pt_ef_andamento_S.fig');
% saveas(fig, '../Figura Simulazioni Articolo/images/pt_ef_andamento_S.fig');