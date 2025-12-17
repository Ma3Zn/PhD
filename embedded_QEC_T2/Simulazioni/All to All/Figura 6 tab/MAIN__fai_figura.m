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

fig_5 = open('images/pt_ef_andamento_S.fig');
pause(0.1);
axx_5 = gca;

fig_6 = open('images/pt_ef_C-phi.fig');
axx_6 = gca;
pause(0.1);

fig_4 = open('images/pt_ef_ripetizione_EC.fig');
pause(0.1);
axx_4 = gca;

fig_4_zoom = open('images/pt_ef_ripetizione_EC_zoom.fig');
pause(0.1);
axx_4_zoom = gca;

% % % fig_4 = open('images/pt_ef_ripetizione_EC_inset.fig');
% % % axx_4_m = fig_4.Children(2);
% % % axx_4_i = fig_4.Children(1); 
% % % pause(0.1);



%% Creiamo la nuova figura unendo tutte le vecchie

figure;
tclf = gcf;
tcl=tiledlayout(2,3);
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

axx_4.Parent=tcl;
axx_4.Layout.Tile=4;
pause(0.1);

inset_size = 0.1;
h_inset = copyobj(axx_4_zoom, tclf);
set(h_inset,'Position', [0.232 0.17 .7*inset_size inset_size])
set(h_inset,'xlim', [2e-5,8e-5], 'ylim',[10e-11,10e-10])
pause(0.1);

axx_5.Parent=tcl;
axx_5.Layout.Tile=5;
pause(0.1);

axx_6.Parent=tcl;
axx_6.Layout.Tile=6;
pause(0.1);

%% Salviamo la figura
saveas(tclf, './images/FINAL_pt_ef_tiled.fig');