% close all
clear
clc

% Script che genera, per varie dimensioni del qudit e vari valori del campo
% esterno da misurare, i grafici delle stime ottenute per Bx e per il
% semiperiodo stimato

% Numero spin da testare
n_spin = 1;

% Valori di spin da testare
spin_vals = 5/2; 

% Numero campi da misurare
n_Bx = 20;

Bx_vals(1,:) = [logspace(-5,-6.70,n_Bx/2) logspace(-6.70, -7.60, 3/2*n_Bx)];

% Per stimare la derivata dT/dB
delta   = 1e-10;
Bx_sens = Bx_vals - Bx_vals*delta;

% Tempi rotazioni da utilizzare
T_rot_vals = 300;

% T2 del sistema [ ns ]
T2 = 5e4;
T1s = [1e6 1e7 1e8 inf];

% livello popolazione da raggiungere per stop ricerca
% lv_pop = 0.7;
lv_pop = 0.5;

% Inizializzazione variabili risultati
T_es        = zeros(n_spin, 2*n_Bx);
T_mis       = zeros(n_spin, 2*n_Bx);
derivata    = zeros(n_spin, 2*n_Bx);
sensitivity = zeros(n_spin, 2*n_Bx);

for T1 = T1s
    T1
    for i = 1:n_spin
    
    %   Spin del sistema da utilizzare
        S = spin_vals(i)
    
    %   Inizializzazione variabili per stima T_max
        T_max = inf;
        Bx_old = Bx_vals(i,1);
    
        for j = 1:2*n_Bx
    
    %%      Primo campo per grafico
    
            Bx = Bx_vals(i,j)
    
    %       Settiamo tempo evoluzione logica tra due istanti di QEC
            T_rot = T_rot_vals(i);
    
    %       Generiamo tutti i dati necessari alla simulazione
            sim = genera_parametri_simulazione(S, Bx, T1, T2, T_rot);
    
    %       Eseguiamo le simulazioni Monte Carlo del protocollo
            tic
            [T_mis(i,j), T_es(i,j)] = esegui_protocollo_sensing(sim, lv_pop, T_max);
    
    %       Controlliamo se il protocollo è riuscito o meno a trovare
    %       l'intersezione
            if (~(T_mis(i,j) > 0))
                break;
            end
    
    %       Pulizia Output
            fprintf('\n');
            toc
    
% % %     %%      Secondo campo per sensitivity
% % %     
% % %     %       Eseguiamo la misura di nu Bx vicino per poter approssimare la
% % %     %       sensitivity
% % %             Bx = Bx_sens(i,j)
% % %     
% % %     %       Generiamo tutti i dati necessari alla simulazione
% % %             sim = genera_parametri_simulazione(S, Bx, T1, T2, T_rot);
% % %     
% % %     %       Eseguiamo le simulazioni Monte Carlo del protocollo
% % %             tic
% % %             T_tmp = esegui_protocollo_sensing(sim, lv_pop, T_max);
     
    %       Settiamo il valore di T_max
            if (j > 1 && j < 2 * n_Bx)
                T_max = T_mis(i,j) * (Bx_vals(i,j) / Bx_vals(i,j+1));
            end
     
% % %     %       Pulizia Output
% % %             fprintf('\n');
% % %             toc
% % %     
% % %     %%      Approssimazione sensitivity
% % %     
% % %     %       Approssimazione derivata
% % %             dSdB = abs( (T_mis(i,j) - T_tmp) / (Bx_vals(i,j) - Bx_sens(i,j)) );
% % %             derivata(i,j) = 1/dSdB;
% % %             sensitivity(i,j) = derivata(i,j) / sqrt(floor(1/(T_es(i,j)*1e-9)));
        end
    end
    
    %%
    
    %   __TTOT__
    %   Salviamo gli esiti della simulazione
    % save(strcat('./sim/__TTOT__', num2str(T_rot_vals(1)), '__', string(T2), '__', string(datetime), '__.mat'));
    % save ris
    
    %%  Plot misure
    
    %   Eseguiamo il plot dei periodi misurati in funzione dei campi effettivi
    fai_plot_mod(T_mis * 1e-9, Bx_vals, spin_vals, false, {'$\tau$ [ s ]', '$B_{x}$ [ T ]', '_Tm_','_Bxe_',num2str(T_rot_vals(1)),num2str(lv_pop),num2str(T_rot_vals(1)),num2str(T1)});
    
    %%  Plot Sensitivity
    
%     fai_plot(Bx_vals, sensitivity, spin_vals(1:n_spin), false, {'$B_x$ [ T ]', '$\eta$ [ THz$^{-1/2}$ ]', '_sensitivity','',num2str(T_rot_vals(1)),num2str(lv_pop),num2str(T_rot_vals(1))});
%     fai_plot(Bx_vals, derivata, spin_vals(1:n_spin), false, {'$B_{x}$', '$|\frac{\partial T}{\partial B_x}|^{-1}$', '_derivata','',num2str(T_rot_vals(1)),num2str(lv_pop),num2str(T_rot_vals(1))});
end