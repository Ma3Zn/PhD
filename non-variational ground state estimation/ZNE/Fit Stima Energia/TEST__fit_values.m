close all
clear
clc

data = "Data/2025-03-20";

backend_noise = "/Brisbane";
backend_noise = "/Fez";
% backend_noise = "/Future";

dim_SbE = "/H_pow_1";

dir = strcat(data, backend_noise, dim_SbE);

load(strcat(dir, "/exp_values_3_tstep.txt"));
load(strcat(dir, "/exp_values_4_tstep.txt"));
load(strcat(dir, "/exp_values_5_tstep.txt"));
load(strcat(dir, "/exp_values_6_tstep.txt"));
load(strcat(dir, "/exp_values_7_tstep.txt"));
load(strcat(dir, "/folding_factors.txt"));
load(strcat(dir, "/exact_exp_val.txt"));

data_range = 1:3:17;
data_range = [1:3:9, 10:17];

exp_values_3_tstep = exp_values_3_tstep(:, data_range);
exp_values_4_tstep = exp_values_4_tstep(:, data_range);
exp_values_5_tstep = exp_values_5_tstep(:, data_range);
exp_values_6_tstep = exp_values_6_tstep(:, data_range);
exp_values_7_tstep = exp_values_7_tstep(:, data_range);
folding_factors = folding_factors(data_range);


% Aggiungiamo una varianza ai dati
prec_est = 0.01;
N_sample = 1e3;
folding_factors = aumenta_ff(folding_factors, N_sample);
exp_values_3_tstep = aggiungi_varianza(exp_values_3_tstep, prec_est, N_sample);
exp_values_4_tstep = aggiungi_varianza(exp_values_4_tstep, prec_est, N_sample);
exp_values_5_tstep = aggiungi_varianza(exp_values_5_tstep, prec_est, N_sample);
exp_values_6_tstep = aggiungi_varianza(exp_values_6_tstep, prec_est, N_sample);
exp_values_7_tstep = aggiungi_varianza(exp_values_7_tstep, prec_est, N_sample);


% % % % Riordiniamo i dati in modo opportuno (se vogliamo combinare FF che
% % % % si overlappano)
% % % [folding_factors, p] = sort(folding_factors);
% % % exp_values_4_tstep = exp_values_4_tstep(:, p);
% % % exp_values_5_tstep = exp_values_5_tstep(:, p);
% % % exp_values_6_tstep = exp_values_6_tstep(:, p);

[n_pow, ~] = size(exp_values_4_tstep);

%%

for i = 1:n_pow
    close all

%%  FIT lineare dei dati
% % %     fit3 = fit(folding_factors', exp_values_3_tstep(i,:)', 'poly1');
% % %     fit4 = fit(folding_factors', exp_values_4_tstep(i,:)', 'poly1');
% % %     fit5 = fit(folding_factors', exp_values_5_tstep(i,:)', 'poly1');
% % %     fit6 = fit(folding_factors', exp_values_6_tstep(i,:)', 'poly1');
% % %     fit7 = fit(folding_factors', exp_values_7_tstep(i,:)', 'poly1');
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
% % %     model = @(x) fit3.p1 * x + fit3.p2; 
% % %     plot(xx, model(xx), 'r', folding_factors, exp_values_3_tstep(i,:), 'b.', 'MarkerSize', 10);
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
% % %     model = @(x) fit7.p1 * x + fit7.p2; 
% % %     plot(xx, model(xx), 'r', folding_factors, exp_values_7_tstep(i,:), 'b.', 'MarkerSize', 10);
% % % 
% % %     plot([0 0 0 0 0], exact_exp_val(i,:) , 'b.', 'MarkerSize', 10);

%%  FIT polinomiale dei dati
% % %     fit4 = fit(folding_factors', exp_values_4_tstep(i,:)', 'poly2');
% % %     fit5 = fit(folding_factors', exp_values_5_tstep(i,:)', 'poly2');
% % %     fit6 = fit(folding_factors', exp_values_6_tstep(i,:)', 'poly2');
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
% % %     model = @(x) fit4.p1 * x.^2 + fit4.p2 * x + fit4.p3; 
% % %     plot(xx, model(xx), 'r', folding_factors, exp_values_4_tstep(i,:), 'b.', 'MarkerSize', 10);
% % % 
% % %     model = @(x) fit5.p1 * x.^2 + fit5.p2 * x + fit5.p3; 
% % %     plot(xx, model(xx), 'r', folding_factors, exp_values_5_tstep(i,:), 'b.', 'MarkerSize', 10);
% % % 
% % %     model = @(x) fit6.p1 * x.^2 + fit6.p2 * x + fit6.p3;
% % %     plot(xx, model(xx), 'r', folding_factors, exp_values_6_tstep(i,:), 'b.', 'MarkerSize', 10);
% % % 
% % %     plot([0 0 0], exact_exp_val(i,:) , 'b.', 'MarkerSize', 10);

%%  FIT esponenziale dei dati
    fit3 = fit(folding_factors', exp_values_6_tstep(i,:)', 'exp1')
    fit4 = fit(folding_factors', exp_values_4_tstep(i,:)', 'exp1')
    fit5 = fit(folding_factors', exp_values_5_tstep(i,:)', 'exp1')
    fit6 = fit(folding_factors', exp_values_6_tstep(i,:)', 'exp1')
    fit7 = fit(folding_factors', exp_values_6_tstep(i,:)', 'exp1')

%%  PLOT FIT
    figure()
    hold on
    grid on

    xx = 0:1e-2:folding_factors(end);

    model = @(x) fit3.a * exp(fit3.b*x); 
    plot(xx, model(xx), folding_factors, exp_values_4_tstep(i,:), 'b.', 'MarkerSize', 10);

    model = @(x) fit4.a * exp(fit4.b*x); 
    plot(xx, model(xx), folding_factors, exp_values_4_tstep(i,:), 'b.', 'MarkerSize', 10);

    model = @(x) fit5.a * exp(fit5.b*x); 
    plot(xx, model(xx), folding_factors, exp_values_5_tstep(i,:), 'b.', 'MarkerSize', 10);

    model = @(x) fit6.a * exp(fit6.b*x); 
    plot(xx, model(xx), folding_factors, exp_values_6_tstep(i,:), 'b.', 'MarkerSize', 10);

    model = @(x) fit7.a * exp(fit7.b*x); 
    plot(xx, model(xx), folding_factors, exp_values_4_tstep(i,:), 'b.', 'MarkerSize', 10);

    plot([0 0 0 0 0], exact_exp_val(i,:) , 'b.', 'MarkerSize', 10);


%%  FIT esponenziale dei dati

% % %     fit3 = fit(folding_factors', exp_values_3_tstep(i,:)', 'exp2')
% % %     fit4 = fit(folding_factors', exp_values_4_tstep(i,:)', 'exp2')
% % %     fit5 = fit(folding_factors', exp_values_5_tstep(i,:)', 'exp2')
% % %     fit6 = fit(folding_factors', exp_values_6_tstep(i,:)', 'exp2')
% % %     fit7 = fit(folding_factors', exp_values_7_tstep(i,:)', 'exp2')
% % % 
% % % %%  PLOT FIT
% % %     figure()
% % %     hold on
% % %     grid on
% % % 
% % %     xx = 0:1e-2:folding_factors(end);
% % % 
% % %     model = @(x) fit3.a * exp(fit3.b*x) + fit3.c * exp(fit3.d*x); 
% % %     plot(xx, model(xx), folding_factors, exp_values_3_tstep(i,:), 'b.', 'MarkerSize', 10);
% % % 
% % %     model = @(x) fit4.a * exp(fit4.b*x) + fit4.c * exp(fit4.d*x); 
% % %     plot(xx, model(xx), folding_factors, exp_values_4_tstep(i,:), 'b.', 'MarkerSize', 10);
% % % 
% % %     model = @(x) fit5.a * exp(fit5.b*x) + fit5.c * exp(fit5.d*x); 
% % %     plot(xx, model(xx), folding_factors, exp_values_5_tstep(i,:), 'b.', 'MarkerSize', 10);
% % % 
% % %     model = @(x) fit6.a * exp(fit6.b*x) + fit6.c * exp(fit6.d*x); 
% % %     plot(xx, model(xx), folding_factors, exp_values_6_tstep(i,:), 'b.', 'MarkerSize', 10);
% % %     
% % %     model = @(x) fit7.a * exp(fit7.b*x) + fit7.c * exp(fit7.d*x); 
% % %     plot(xx, model(xx), folding_factors, exp_values_7_tstep(i,:), 'b.', 'MarkerSize', 10);
% % % 
% % %     plot([0 0 0 0 0], exact_exp_val(i,:) , 'b.', 'MarkerSize', 10);

    pause()
end