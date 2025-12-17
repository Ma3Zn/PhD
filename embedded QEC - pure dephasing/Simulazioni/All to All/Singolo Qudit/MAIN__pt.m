% Script per il plot delle pseudo-trhesholds per UN'operazione logica
% partendo da UN unico stato iniziale per un sistema con connessione tutti
% con tutti

%% Pulizia ambiente di lavoro
clc

warning('off','MATLAB:legend:IgnoringExtraEntries');
warning('off','MATLAB:plot:IgnoreImaginaryXYPart');
warning('off','MATLAB:logm:nonPosRealEig');

%% Settiamo i parametri del grafico

% Apriamo la figura su cui fare il grafico
fig = figure(1);

% Tempo IDLE da simulare (in ns)
T_IDLE = 10;
% T_IDLE = 10 * T_IDLE;

% Numero di valori di T2 da testare
N_T2 = 50;

% Range di valori di T2 da testare (in ns)
T2_vals = logspace(0,8,N_T2);
T2_vals = T2_vals(end:-1:1);

% T2_vals = 5e4;

% Range di spin di cui fare il grafico
dim_vals = 4:2:6;

%% Facciamo il grafico della soluzione benchmark (spin 1/2 non corretto)

grafico_errore_1_2;

%% Facciamo i grafici della procedura di EC per tutti gli spin da testare

for dim = dim_vals
    tic;
    dim

    [errore_medio, T] = calcola_errore_medio(dim, T2_vals);

    fai_grafico(errore_medio, T2_vals, 0);
    toc;
end

% Salviamo il grafico
legend('S = 1/2' , 'lv = 4' , 'lv = 6' , 'lv = 8' , 'lv = 10' , 'lv = 12', ...
       'Location', 'southeast');
xlabel('$1/T_2$','Interpreter','latex','FontSize',18);
ylabel('$\mathcal{E} = 1 - \mathcal{F}^2$','Interpreter','latex','FontSize',18);
saveas(fig, './images/pt.fig');