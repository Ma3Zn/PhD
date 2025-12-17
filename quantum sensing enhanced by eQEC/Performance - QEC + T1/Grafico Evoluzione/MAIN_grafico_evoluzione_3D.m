% close all
clear
clc

% Script che genera, per varie dimensioni del qudit e vari valori del campo
% esterno da misurare, i grafici dell'evoluzione delle popolazioni logiche

% Numero spin da testare
n_spin = 3;

% Valori di spin da testare
spin_vals = [3/2 5/2 7/2];
% spin_vals = [7/2];

% Stabiliamo se vogliamo o meno una error correction perfetta
perfect_EC = true;

% Valori del campo da misurare
Bx_vals = cell(1,3);
Bx_vals{1} = [5e-7 1e-7 5e-8];
Bx_vals{2} = [5e-7 1e-7 5e-8];
Bx_vals{3} = [5e-7 1e-7 5e-8];
% Bx_vals{1} = [5e-6 1e-6 5e-7 1e-7 5e-8];
% Bx_vals{1} = [1e-4 1e-5 1e-6];

% Tempi rotazioni da utilizzare
T_rot_vals = [400 400 400];

% T2 del sistema [ ns ]
T2 = 5e4;

%   Generazione figura per rispettivo spin
figure()
grid on
hold on

for i = 1:n_spin

%   Spin del sistema da utilizzare
    S = spin_vals(i)

    for j = 1:length(Bx_vals{i})

%%      Primo campo per grafico

        Bx = Bx_vals{i}(j)

%       Settiamo tempo evoluzione logica tra due istanti di QEC
        T_rot = T_rot_vals(i);

%       Generiamo tutti i dati necessari alla simulazione
        sim = genera_parametri_simulazione(S, Bx, T2, T_rot);

%       Eseguiamo le simulazioni Monte Carlo del protocollo
        tic
        sim = esegui_protocollo_sensing(sim);

%       Pulizia Output
        fprintf('\n');
        toc

%%      Plot evoluzione
        plot3(linspace(0,sim.T_fin,sim.it-1) * 1e-9, S*ones(length(sim.pop),1), sim.pop, 'LineWidth', 2.5);
    end
end

view(3)