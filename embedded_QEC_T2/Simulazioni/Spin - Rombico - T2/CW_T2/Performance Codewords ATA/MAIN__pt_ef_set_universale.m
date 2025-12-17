% Script per il plot della pseudo-trheshold mediata per un SET UNIVERSALE
% di operazioni logiche a singolo corpo ottenute tramite la computazione
% dell'entaglement fidelity

%% Pulizia ambiente di lavoro
clear
clc

warning('off','MATLAB:legend:IgnoringExtraEntries');
warning('off','MATLAB:plot:IgnoreImaginaryXYPart');
warning('off','MATLAB:logm:nonPosRealEig');

%% Settiamo i parametri del grafico

% Apriamo la figura su cui fare il grafico
fig = figure(1);

% Tempo gate fisico da simulare (in ns)
T_gate_fisico = 10;

% Numero di valori di T2 da testare
N_T2 = 20;

% Range di valori di T2 da testare (in ns)
T2_vals = logspace(0,8,N_T2);
T2_vals = T2_vals(end:-1:1);
% T2_vals = 1e7;

% Range dimensioni di cui fare il grafico
dim_vals = 4:2:10;
dim_vals = 8;

% Valori angoli per un set di rotazioni piane universali (singolo corpo)
angle_vals = [pi/4 pi; pi/2 pi; pi/2 -pi/2; pi/2 -pi/4; pi/2 -pi/8]';
angle_vals = [0 0]';

% Possibili stati iniziali del qubit logico
init_superposition = [1 0; 0 1; 1/sqrt(2) 1/sqrt(2); 1/sqrt(2) -1/sqrt(2); 1/sqrt(2) 1i/sqrt(2); 1/sqrt(2) -1i/sqrt(2)]';
init_superposition = [1/sqrt(2) -1i/sqrt(2)]';

%% Facciamo il garfico di un qubit NON corretto

dev = grafico_errore_medio_1_2(T2_vals, angle_vals, init_superposition);

%% Facciamo il grafico di un qubit logico

for dim = dim_vals
    tic;
%     figure();

    dim
    tmp = grafico_errore_medio(dim, T2_vals, angle_vals, init_superposition);

%   Concateniamo le deviazioni della media
    dev = [dev; tmp];

    toc;

    % Mettiamo legenda e label
    title(strcat(num2str(dim),'lv'));
    legend('lv$ = 2$','lv $ = 4$', 'lv $ = 6$', 'lv $ = 8$', 'lv $ = 10$', 'lv $ = 12$', 'lv $ = 14$', 'lv $ = 16$', 'Interpreter','latex','Location', 'southeast','FontSize',16.2);
    
    xlabel('$1/T_2$','Interpreter','latex','FontSize',18);
    ylabel('$\mathcal{E} = \overline{1 - \mathcal{F}_{e}{}^2}$','Interpreter','latex','FontSize',18);

    % Salviamo il grafico
    saveas(fig, strcat('./images/pt_ef_',num2str(dim),'lv.fig'));

end

% plottiamo gli intorni tubolari
% plot(1./T2_vals, pulisci_matrice(dev, 1e-2 * eps), '--k', 'HandleVisibility','off');