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
% T2_vals = [1e5, 1e4, 1e3];
T2_vals = 5e4

% Valori angoli per un set di rotazioni piane universali (singolo corpo)
% angle_vals = [pi/4 pi; pi/2 pi; pi/2 -pi/2; pi/2 -pi/4; pi/2 -pi/8]';
angle_vals = [0 0]';

% Possibili stati iniziali del qubit logico
init_superposition = [1 0; 0 1; 1/sqrt(2) 1/sqrt(2); 1/sqrt(2) -1/sqrt(2); 1/sqrt(2) 1i/sqrt(2); 1/sqrt(2) -1i/sqrt(2)]';
% init_superposition = [1/sqrt(2) -1i/sqrt(2)]';

% Range di spin di cui fare il grafico
dim_vals = 4:2:6;

% Iterazioni da testare
iter = 1:6;
% iter = 1:4;

%% per generare correttamente l'input del programma
genera_input_1_2;

%% Facciamo i grafici della procedura di EC per tutti gli spin da testare
errore_medio = zeros(length(iter),1);

for T2 = T2_vals
    T2
    
    tic
    for dim = dim_vals
        dim

        idx1 = 1;
        for it = iter
            it
            errore_medio(idx1) = grafico_errore_medio(dim, T2, angle_vals, init_superposition, it);
            idx1 = idx1 + 1;
        end

        errore_medio = pulisci_matrice(errore_medio, 1e3*eps);
        plot(iter, errore_medio, '-x', 'MarkerSize', 10, 'LineWidth', 3);
        grid on
        hold on
    end
    toc
end

%%

% Salviamo il grafico
legend('lv = 4', 'lv = 6', 'lv = 8', 'Location', 'southeast');
xlabel('\textit{lv}','Interpreter','latex','FontSize',18);
ylabel('$\mathcal{E} = 1 - \mathcal{F}^2$','Interpreter','latex','FontSize',18);

saveas(fig, './images/pt_ef_andamento_ripetizioni_EC.fig');
saveas(fig, '../Figura 4 tab MAIN/images/pt_ef_andamento_ripetizioni_EC.fig');