close all
clear
clc

data = "Data/28-2-25";

load(strcat(data, "/exp_values_4_tstep.txt"));
load(strcat(data, "/exp_values_5_tstep.txt"));
load(strcat(data, "/exp_values_6_tstep.txt"));
load(strcat(data, "/folding_factors.txt"));
load(strcat(data, "/exact_exp_val.txt"));

data_range = 1:9;
% data_range = [1:5, 10:2:18];

exp_values_4_tstep = exp_values_4_tstep(:, data_range);
exp_values_5_tstep = exp_values_5_tstep(:, data_range);
exp_values_6_tstep = exp_values_6_tstep(:, data_range);
folding_factors = folding_factors(data_range);

% Riordiniamo i dati in modo opportuno
[folding_factors, p] = sort(folding_factors);
exp_values_4_tstep = exp_values_4_tstep(:, p);
exp_values_5_tstep = exp_values_5_tstep(:, p);
exp_values_6_tstep = exp_values_6_tstep(:, p);

[n_pow, ~] = size(exp_values_4_tstep);

%%

for i = 1:n_pow
    close all

% % % %%  FIT lineare dei dati
% % %     fit4 = fit(folding_factors', exp_values_4_tstep(i,:)', 'poly1');
% % %     fit5 = fit(folding_factors', exp_values_5_tstep(i,:)', 'poly1');
% % %     fit6 = fit(folding_factors', exp_values_6_tstep(i,:)', 'poly1');
% % % 
% % % %%  PLOT FIT
% % %     figure()
% % %     hold on
% % %     grid on
% % % 
% % % %     plot(fit4, folding_factors, exp_values_4_tstep(i,:));
% % % %     plot(fit5, folding_factors, exp_values_5_tstep(i,:));
% % % %     plot(fit6, folding_factors, exp_values_6_tstep(i,:));
% % % 
% % %     xx = 0:1e-2:folding_factors(end);
% % % 
% % %     model = @(x) fit4.p1 * x + fit4.p2;
% % %     plot(xx, model(xx), 'r', folding_factors, exp_values_4_tstep(i,:), 'b.', 'MarkerSize', 10);
% % % 
% % %     model = @(x) fit5.p1 * x + fit5.p2;
% % %     plot(xx, model(xx), 'r', folding_factors, exp_values_5_tstep(i,:), 'b.', 'MarkerSize', 10);
% % % 
% % %     model = @(x) fit6.p1 * x + fit6.p2;
% % %     plot(xx, model(xx), 'r', folding_factors, exp_values_6_tstep(i,:), 'b.', 'MarkerSize', 10);
% % % 
% % %     plot([0 0 0], exact_exp_val(i,:) , 'b.', 'MarkerSize', 10);

%%  FIT polinomiale dei dati
    fit4 = fit(folding_factors', exp_values_4_tstep(i,:)', 'poly2');
    fit5 = fit(folding_factors', exp_values_5_tstep(i,:)', 'poly2');
    fit6 = fit(folding_factors', exp_values_6_tstep(i,:)', 'poly2');

%%  PLOT FIT
    figure()
    hold on
    grid on

%     plot(fit4, folding_factors, exp_values_4_tstep(i,:));
%     plot(fit5, folding_factors, exp_values_5_tstep(i,:));
%     plot(fit6, folding_factors, exp_values_6_tstep(i,:));

    xx = 0:1e-2:folding_factors(end);

    model = @(x) fit4.p1 * x.^2 + fit4.p2 * x + fit4.p3; 
    plot(xx, model(xx), 'b', folding_factors, exp_values_4_tstep(i,:), 'bo', 'MarkerSize', 10, LineWidth=2.5);

    model = @(x) fit5.p1 * x.^2 + fit5.p2 * x + fit5.p3; 
    plot(xx, model(xx), 'r', folding_factors, exp_values_5_tstep(i,:), 'ro', 'MarkerSize', 10, LineWidth=2.5);

    model = @(x) fit6.p1 * x.^2 + fit6.p2 * x + fit6.p3;
    plot(xx, model(xx), 'k', folding_factors, exp_values_6_tstep(i,:), 'ko', 'MarkerSize', 10, LineWidth=2.5);

    plot(0, exact_exp_val(i,1) , 'bo', 'MarkerSize', 10, LineWidth=2.5);
    plot(0, exact_exp_val(i,2) , 'ro', 'MarkerSize', 10, LineWidth=2.5);
    plot(0, exact_exp_val(i,3) , 'ko', 'MarkerSize', 10, LineWidth=2.5);

    set(gca,'FontSize',15);
    title('ZNE of $\langle H \rangle$ - Quadratic Model', 'Interpreter', 'latex', 'FontSize',20);
    xlabel('Gate Folding Factor','Interpreter','latex',"FontSize",20); 
    ylabel('E','Interpreter','latex',"FontSize",20,'Rotation',0);

%%  FIT esponenziale dei dati
    fit4 = fit(folding_factors', exp_values_4_tstep(i,:)', 'exp1');
    fit5 = fit(folding_factors', exp_values_5_tstep(i,:)', 'exp1');
    fit6 = fit(folding_factors', exp_values_6_tstep(i,:)', 'exp1');

%%  PLOT FIT
    figure()
    hold on
    grid on

    xx = 0:1e-2:folding_factors(end);

    model = @(x) fit4.a * exp(fit4.b*x); 
    plot(xx, model(xx), 'b', folding_factors, exp_values_4_tstep(i,:), 'bo', 'MarkerSize', 10, LineWidth=2.5);

    model = @(x) fit5.a * exp(fit5.b*x); 
    plot(xx, model(xx), 'r', folding_factors, exp_values_5_tstep(i,:), 'ro', 'MarkerSize', 10, LineWidth=2.5);

    model = @(x) fit6.a * exp(fit6.b*x); 
    plot(xx, model(xx), 'k', folding_factors, exp_values_6_tstep(i,:), 'ko', 'MarkerSize', 10, LineWidth=2.5);

    plot(0, exact_exp_val(i,1) , 'bo', 'MarkerSize', 10, LineWidth=2.5);
    plot(0, exact_exp_val(i,2) , 'ro', 'MarkerSize', 10, LineWidth=2.5);
    plot(0, exact_exp_val(i,3) , 'ko', 'MarkerSize', 10, LineWidth=2.5);

    set(gca,'FontSize',15);
    title('ZNE of $\langle H \rangle$ - Exponential Model', 'Interpreter', 'latex', 'FontSize',20);
    xlabel('Gate Folding Factor','Interpreter','latex',"FontSize",20);

%%  FIT esponenziale dei dati
    fit4 = fit(folding_factors', exp_values_4_tstep(i,:)', 'exp2')
    fit5 = fit(folding_factors', exp_values_5_tstep(i,:)', 'exp2')
    fit6 = fit(folding_factors', exp_values_6_tstep(i,:)', 'exp2')

%%  PLOT FIT
    figure()
    hold on
    grid on

    xx = 0:1e-2:folding_factors(end);

    model = @(x) fit4.a * exp(fit4.b*x) + fit4.c * exp(fit4.d*x); 
    plot(xx, model(xx), 'b', folding_factors, exp_values_4_tstep(i,:), 'bo', 'MarkerSize', 10, LineWidth=2.5);

    model = @(x) fit5.a * exp(fit5.b*x) + fit5.c * exp(fit5.d*x); 
    plot(xx, model(xx), 'r', folding_factors, exp_values_5_tstep(i,:), 'ro', 'MarkerSize', 10, LineWidth=2.5);

    model = @(x) fit6.a * exp(fit6.b*x) + fit6.c * exp(fit6.d*x); 
    plot(xx, model(xx), 'k', folding_factors, exp_values_6_tstep(i,:), 'ko', 'MarkerSize', 10, LineWidth=2.5);

    plot(0, exact_exp_val(i,1) , 'bo', 'MarkerSize', 10, LineWidth=2.5);
    plot(0, exact_exp_val(i,2) , 'ro', 'MarkerSize', 10, LineWidth=2.5);
    plot(0, exact_exp_val(i,3) , 'ko', 'MarkerSize', 10, LineWidth=2.5);

    set(gca,'FontSize',15);
    title('ZNE of $\langle H \rangle$ - Bi-Exponential Model', 'Interpreter', 'latex', 'FontSize',20);
    xlabel('Gate Folding Factor','Interpreter','latex',"FontSize",20);

    pause()
end