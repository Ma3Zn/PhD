close all
clc

%% FIG 2-(a)

% Recuperiamo i dati del plot
fig = openfig("./Singolo Qudit/GS_pt_ef.fig");

axObjs = fig.Children;
dataObjs = axObjs(2).Children;

x = dataObjs(1).XData;

y1 = dataObjs(1).YData;
y2 = dataObjs(2).YData;
y3 = dataObjs(3).YData;
y4 = dataObjs(4).YData;
y5 = dataObjs(5).YData;
y6 = dataObjs(6).YData;

data_fig2_a = [x; y1; y2; y3; y4; y5; y6].'

%% FIG 2-(a) -- Inset
% Recuperiamo i dati del plot
fig = openfig("./Singolo Qudit/BP_pt_ef.fig");

axObjs = fig.Children;
dataObjs = axObjs(2).Children;

x = dataObjs(1).XData;

y1 = dataObjs(1).YData;
y2 = dataObjs(2).YData;
y3 = dataObjs(3).YData;
y4 = dataObjs(4).YData;
y5 = dataObjs(5).YData;
y6 = dataObjs(6).YData;

data_fig2_a__inset = [x; y1; y2; y3; y4; y5; y6].'

%% FIG 2-(b)

% Recuperiamo i dati del plot
fig = openfig("./Singolo Qudit/GS_pt_ef_andamento_S.fig");

axObjs = fig.Children;
dataObjs = axObjs(2).Children;

x = dataObjs(1).XData;

y1 = dataObjs(1).YData;
y2 = dataObjs(2).YData;
y3 = dataObjs(3).YData;

% Recuperiamo i dati del plot
fig = openfig("./Singolo Qudit/BP_pt_ef_andamento_S.fig");

axObjs = fig.Children;
dataObjs = axObjs(2).Children;

y4 = dataObjs(1).YData;
y5 = dataObjs(2).YData;
y6 = dataObjs(3).YData;

data_fig2_b = [x; y1; y2; y3; y4; y5; y6].'

%% Fig 3
% Recuperiamo i dati del plot
fig = openfig("./Doppio Qudit/GS_pt.fig");

axObjs = fig.Children;
dataObjs = axObjs(2).Children;

x1 = dataObjs(1).XData;
y1 = dataObjs(1).YData;
y2 = dataObjs(2).YData;
y3 = dataObjs(3).YData;

% Recuperiamo i dati del plot
fig = openfig("./Doppio Qudit/BP_pt.fig");

axObjs = fig.Children;
dataObjs = axObjs(2).Children;

y4 = dataObjs(1).YData;

x2 = dataObjs(2).XData;
y5 = dataObjs(2).YData;

data_fig3 = [x1; y1; y2; y3; y4].'
data_fig3_4 = [x2; y5].'

%% FIG mis 4
fig = openfig("./Errore Misura/errore_misura_rep_4_FINAL.fig");
axObjs = fig.Children;
dataObjs = axObjs(2).Children;

x = dataObjs(1).XData;

y1 = dataObjs(1).YData;
y2 = dataObjs(2).YData;
y3 = dataObjs(3).YData;
y4 = dataObjs(4).YData;
y5 = dataObjs(5).YData;

data_fig_mis_4 = [x; y5; y1; y2; y3; y4].'

%% FIG mis 6
fig = openfig("./Errore Misura/errore_misura_rep_6_FINAL.fig");
axObjs = fig.Children;
dataObjs = axObjs(2).Children;

x = dataObjs(1).XData;

y1 = dataObjs(1).YData;
y2 = dataObjs(2).YData;
y3 = dataObjs(3).YData;
y4 = dataObjs(4).YData;
y5 = dataObjs(5).YData;

data_fig_mis_6 = [x; y5; y3; y2; y1; y4].'

%% Save

save ./DATA_SHARE/data_fig2_a.txt data_fig2_a -ascii
save ./DATA_SHARE/data_fig2_a__inset.txt data_fig2_a__inset -ascii
save ./DATA_SHARE/data_fig2_b.txt data_fig2_b -ascii
save ./DATA_SHARE/data_fig3.txt data_fig3 -ascii
save ./DATA_SHARE/data_fig3_4.txt data_fig3_4 -ascii
save ./DATA_SHARE/data_fig_mis_4.txt data_fig_mis_4 -ascii
save ./DATA_SHARE/data_fig_mis_6.txt data_fig_mis_6 -ascii