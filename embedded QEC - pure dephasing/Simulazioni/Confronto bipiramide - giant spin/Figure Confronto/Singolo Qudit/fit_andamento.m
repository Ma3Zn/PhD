close all
clc

% Recuperiamo i dati del plot
fig = openfig("GS_pt_ef_andamento_S.fig");

axObjs = fig.Children;
dataObjs = axObjs(2).Children;

x1 = dataObjs(1).XData;
y1 = dataObjs(1).YData;

x2 = dataObjs(2).XData;
y2 = dataObjs(2).YData;

x3 = dataObjs(3).XData;
y3 = dataObjs(3).YData;

% Rimuoviamo i punti che non ci garbano
x1 = x1(2:end);
y1 = y1(2:end);

x2 = x2(2:end);
y2 = y2(2:end);

x3 = x3(2:end-2);
y3 = y3(2:end-2);

p1 = polyfit(x1, log(y1), 1);
p2 = polyfit(x2, log(y2), 1);
p3 = polyfit(x3, log(y3), 1);

% Proviamo a fare qualche tipo di fit
xx = linspace(2,12,1e3);

semilogy(xx, exp(polyval(p1, xx)) ,'LineWidth', 3);
semilogy(xx, exp(polyval(p2, xx)) ,'LineWidth', 3);
semilogy(xx, exp(polyval(p3, xx)) ,'LineWidth', 3);

% % % p1 = polyfit(x1, y1, 2);
% % % p2 = polyfit(x2, y2, 2);
% % % p3 = polyfit(x3, y3, 2);
% % % 
% % % semilogy(xx, polyval(p1, xx) ,'LineWidth', 3);
% % % semilogy(xx, polyval(p2, xx) ,'LineWidth', 3);
% % % semilogy(xx, polyval(p3, xx) ,'LineWidth', 3);

% % % % % % % % % % % % % % semilogy(xx, polyval([-1,zeros(1,1e1-1)],xx));
% % % % % % % % % % % % % % semilogy(xx, polyval([-1,zeros(1,1e2-1)],xx));
% % % % % % % % % % % % % % semilogy(xx, polyval([-1,zeros(1,1e3-1)],xx));
