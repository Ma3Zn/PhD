close all
clear
clc

%% Recuperiamo le figure
fig_1 = open('images/pt_ef_errori_totali.fig');
pause(0.1);
axx_1 = gca;

fig_2 = open('images/pt_ef_errori_NON_correggibili_singolo.fig');
pause(0.1);
axx_2 = gca;

fig_3 = open('images/pt_ef_tomografia_logica.fig');
pause(0.1);
axx_3 = gca;





%% Creiamo la nuova figura unendo tutte le vecchie

figure;
tclf = gcf;
tcl=tiledlayout(1,3);
pause(0.1);

axx_1.Parent=tcl;
axx_1.Layout.Tile=1;
pause(0.1);

axx_2.Parent=tcl;
axx_2.Layout.Tile=2;
pause(0.1);

axx_3.Parent=tcl;
axx_3.Layout.Tile=3;
pause(0.1);


%% Salviamo la figura
saveas(tclf, './images/FINAL_pt_ef_tiled.fig');