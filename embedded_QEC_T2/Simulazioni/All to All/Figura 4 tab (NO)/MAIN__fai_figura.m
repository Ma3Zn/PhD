close all
clear
clc

%% Recuperiamo le figure

fig_1 = open('images/pt_ef_Tomografia_logica_allungamento.fig');
pause(0.01);
ax_1 = gca;

fig_2 = open('images/pt_ef_C-phi.fig');
pause(0.01);
ax_2 = gca;

fig_3 = open('images/pt_ef_andamento_S.fig');
pause(0.01);
ax_3 = gca;

fig_4 = open('images/pt_ef_andamento_ripetizioni_EC_4.fig');
pause(0.01);
ax_4 = gca;

fig_5 = open('images/pt_ef_andamento_ripetizioni_EC_6.fig');
pause(0.01);
ax_5 = gca;

%% Creiamo la nuova figura unendo tutte le vecchie

figure;
sbp(1) = subplot(5,4,[1 2 5 6],ax_1);
sbp(2) = subplot(5,4,[3 4 7 8],ax_2);
sbp(3) = subplot(5,4,[9 10],ax_3);
sbp(4) = subplot(5,4,11,ax_4);
sbp(5) = subplot(5,4,12,ax_5);

% Paste figures on the subplots

% Add legends
% l(1)=legend(h(1))
% l(2)=legend(h(2))


%% Salviamo la figura
saveas(sbp, './images/FINAL_pt_ef_tiled.fig');