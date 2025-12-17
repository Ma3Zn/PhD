% close all
clear
clc

% Script che genera, per varie dimensioni del qudit e vari valori del campo
% esterno da misurare, i grafici delle stime ottenute per Bx e per il
% semiperiodo stimato

% Numero spin da testare
n_spin = 2;

% Valori di spin da testare
spin_vals = [3/2 5/2 7/2];
% spin_vals = [7/2];

% Numero iterate algoritmo di ricerca da eseguire
n_search = 8;

% Guess iniziali campi da misurare [T]
Bx_vals(1) = 2.1e-7;
Bx_vals(2) = 3.5e-8;
Bx_vals(3) = 1.5e-8;

% Tempi massimi per i Guess dei tempi iniziali [ns]
T_max_vals(1) = 0.015e9;
T_max_vals(2) = 0.450e9;
T_max_vals(3) = 1.000e9;

% Tempo di rotazione del Guess iniziale
T_rot_guess = 500;

% Numero di tempi di rotazione da utilizzare
n_rot = 20;

% Tempi rotazioni da utilizzare
T_rot_vals = logspace(log10(150),5,n_rot);

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

%       Recuperiamo il tempo di rotazione attuale
        T_rot = T_rot_vals(j)

%       Recuperiamo il Guess scalato opportunamente
%         Bx = Bx_vals(i) * T_rot / (4^log10(T_rot_vals(i)) * T_rot_guess);
        Bx = Bx_vals(i);

%       Stima iniziale T_max
        T_max = T_max_vals(i);

%       Verifica del campo attuale
        tic
        misurabile = verifica_detectabilita_campo(S, Bx, T2, T_rot, lv_pop, T_max);

%       In base ad esito verifica del guess iniziale individuiamo il primo
%       intervallo in cui eseguire la ricerca dicotomica del campo minimo
%       misurabile
        if (misurabile)
%           Cerchiamo l'intervallo in avanti
            [Bx_a, Bx_b, T_max] = cerca_intervallo_avanti(S, Bx, T2, T_rot, lv_pop, T_max);
        else
%           Cerchiamo l'intervallo all'indietro 
%           [ATTENZIONE:: NON SOVRASCRIVERE T_MAX IN QUESTA SITUAZIONE]
            [Bx_a, Bx_b] = cerca_intervallo_indietro(S, Bx, T2, T_rot, lv_pop, T_max);
        end

%       Una volta che abbiamo il nostro intervallo eseguiamo una ricerca
%       dicotomica del più piccolo campo misurabile dal protocollo [TBD]
        Bx_min_est(i,j) = ricerca_dicotomica_Bx_min(S, Bx_a, Bx_b, T2, T_rot, lv_pop, T_max, n_search);
        toc
    end
end

%%

%   __TTOT__
%   Salviamo gli esiti della simulazione
save(strcat('./sim/__BXMIN__', string(T2), '__', string(datetime), '__.mat'));

%%  Plot stime

%   Eseguiamo il plot dei periodi misurati in funzione dei campi effettivi
fai_plot([T_rot_vals;T_rot_vals;T_rot_vals], Bx_min_est, spin_vals(1:n_spin), false, {'Intervallo QEC [ $\mu$s ]', '$B_{x}$ [ T ]', '_BXmin_', '_Trot_', num2str(T_rot_vals(1)), num2str(lv_pop), num2str(T_rot_vals(1)*100)});
