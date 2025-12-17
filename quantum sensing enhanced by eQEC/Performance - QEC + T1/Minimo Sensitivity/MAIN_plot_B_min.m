% close all
clear
clc

% Script che genera, per varie dimensioni del qudit e vari valori del campo
% esterno da misurare, i grafici delle stime ottenute per Bx e per il
% semiperiodo stimato

% Numero spin da testare
n_spin = 2;

% Valori di spin da testare
spin_vals = [3/2 5/2];

% Numero iterate algoritmo di ricerca da eseguire
n_search = 8;

% Guess iniziali campi da misurare [T]
Bx_vals(1) = 1e-5;
Bx_vals(2) = 6e-6;

% Tempi massimi per i Guess dei tempi iniziali [ns]
T_max_vals(1) = 0.01e9;
T_max_vals(2) = 0.07e9;

% Tempo di rotazione del Guess iniziale
T_rot_guess = 500;

% Numero di tempi di rotazione da utilizzare
n_rot = 20;

% Tempi rotazioni da utilizzare
T_rot_vals = logspace(2,4,n_rot);

% T2 del sistema [ ns ]
T2 = 5e4;
T1 = 1e8;

% livello popolazione da raggiungere per stop ricerca
% lv_pop = 0.7;
lv_pop = 0.5;

% Inizializzazione varibaile contenente i campi minimi stimati
Bx_min_est = zeros(n_spin, n_rot);

% Ciclo sui vari spin
for i = 1:n_spin

%   Spin del sistema da utilizzare
    S = spin_vals(i)

%   Inizializzazione variabili per stima T_max
    T_max = inf;
    Bx_old = Bx_vals(i);

%   Ciclo tempi di rotazione
    for j = 1:n_rot
        tic
%       Recuperiamo il tempo di rotazione attuale
        T_rot = T_rot_vals(j)

%       Recuperiamo il valore iniziale del campo
        Bx = Bx_vals(i);

%       Stima iniziale T_max
        T_max = T_max_vals(i);

%       Inizializzazione varibile sensitività
        sens     = 1;

%       Eseguiamo la ricerca in avanti adattiva del minimo
        for k = 1:8
%           Calcoliamo lo step della ricerca
            step = 1/(10^k);

%           In base ad esito verifica del guess iniziale individuiamo il primo
%           intervallo in cui eseguire la ricerca dicotomica del campo minimo
%           misurabile
            [Bx, sens, T_max] = cerca_in_avanti(S, Bx, sens, step, T2, T_rot, lv_pop, T_max);
        end
        toc

%       Salviamo la stima del campo minimo attuale
        Bx_min_est(i,j) = Bx;

%       Salviamo lo stato attuale
        save ris
    end
end

%%

%   __TTOT__
% Salviamo gli esiti della simulazione

% save(strcat('./sim/__BXMIN__', string(T2), '__', string(datetime), '__.mat'));

%%  Plot stime

% Eseguiamo il plot dei periodi misurati in funzione dei campi effettivi
fai_plot([T_rot_vals;T_rot_vals;T_rot_vals], Bx_min_est, spin_vals(1:n_spin), false, {'Intervallo QEC [ $\mu$s ]', '$B_{x}$ [ T ]', '_BXmin_', '_Trot_', num2str(T_rot_vals(1)), num2str(lv_pop), num2str(T_rot_vals(1)*100)});
