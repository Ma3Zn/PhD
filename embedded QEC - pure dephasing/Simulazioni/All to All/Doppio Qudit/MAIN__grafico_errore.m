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
T_gate_fisico = 2;

% Numero di valori di T2 da testare
N_T2 = 20;

% Range di valori di T2 da testare (in ns)
% T2_vals = logspace(0,8,N_T2);
T2_vals = [logspace(0,4,8), 2e5, 1e7];
T2_vals = T2_vals(end:-1:1);
% T2_vals = 5e4;

% Range dimensioni di cui fare il grafico
% dim_vals = 4:2:6;
dim_vals = 4;

% Angoli rotazione
angle_vals = pi;

% Possibili stati iniziali dei qubit logici
init_superposition = [1 0; 0 1; 1/sqrt(2) 1/sqrt(2); 1/sqrt(2) -1/sqrt(2); 1/sqrt(2) 1i/sqrt(2); 1/sqrt(2) -1i/sqrt(2)].';
% init_superposition = [1/sqrt(2) -1i/sqrt(2)]';

%% Facciamo il garfico di un qubit NON corretto

dev = grafico_errore_medio_1_2(T2_vals, angle_vals, init_superposition);

%% Facciamo il grafico di un qubit logico

for dim = dim_vals

    fprintf('INIZIO simulazione %d livelli\n\n', dim);

    tmp = grafico_errore_medio(fig, dim, T2_vals, angle_vals, init_superposition);

    fprintf('FINE simulazione %d livelli\n', dim);

%   Concateniamo le deviazioni della media
    dev = [dev; tmp];
end

% Mettiamo legenda e label
legend('S = 1/2' , 'lv = 4' , 'lv = 6' , 'lv =8' , 'lv = 10' , 'lv = 12', ...
       'Location', 'southeast');
xlabel('$1/T_2$','Interpreter','latex','FontSize',18);
ylabel('$\mathcal{E} = \overline{1 - \mathcal{F}_{e}{}^2}$','Interpreter','latex','FontSize',18);

%%

% plottiamo gli intorni tubolari
plot(1./T2_vals, pulisci_matrice(dev, 1e-2 * eps), '--k', 'HandleVisibility','off');

% Salviamo il grafico
saveas(fig, './images/pseudo_thresholds.fig');
