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
% T2_vals = 5e4;

% Range dimensioni di cui fare il grafico
dim_vals = 4:2:12;
dim_vals = 4;

% Valori angoli per un set di rotazioni piane universali (singolo corpo)
angle_vals = [pi/4 pi; pi/2 pi; pi/2 -pi/2; pi/2 -pi/4; pi/2 -pi/8]';
% angle_vals = [pi/4 pi]';

% Possibili stati iniziali del qubit logico
init_superposition = [1 0; 0 1; 1/sqrt(2) 1/sqrt(2); 1/sqrt(2) -1/sqrt(2); 1/sqrt(2) 1i/sqrt(2); 1/sqrt(2) -1i/sqrt(2)]';
% init_superposition = [1/sqrt(2) -1i/sqrt(2)]';

%% Facciamo il garfico di un qubit NON corretto

dev = grafico_errore_medio_1_2(T2_vals, angle_vals, init_superposition);

%% Facciamo il grafico di un qubit logico

for dim = dim_vals
    tic;
    
    dim

    tmp = grafico_errore_medio(dim, T2_vals, angle_vals, init_superposition);

%   Concateniamo le deviazioni della media
    dev = [dev; tmp];
    
    toc;
end

% Mettiamo legenda e label
% legend('S = 1/2' , 'lv = 4' , 'lv = 6' , 'lv = 8' , 'lv = 10' , 'lv = 12', ...
%        'Location', 'southeast');
xlabel('$1/T_2$','Interpreter','latex','FontSize',18);
ylabel('$\mathcal{E}$','Interpreter','latex','FontSize',18,'Rotation',0);
% ylabel('$\mathcal{E} = 1 - \mathcal{F}_{e}{}^2$','Interpreter','latex','FontSize',18);

% plottiamo gli intorni tubolari
% plot(1./T2_vals, pulisci_matrice(dev, 1e-2 * eps), '--k', 'HandleVisibility','off');

%%
% Salviamo il grafico
saveas(fig, './images/pt_ef.fig');

%%
% saveas(fig, './images/pt_ef_Tomografia_logica_allungamento.fig');
% saveas(fig, '../Figura 4 tab MAIN/images/pt_ef_Tomografia_logica_allungamento.fig');