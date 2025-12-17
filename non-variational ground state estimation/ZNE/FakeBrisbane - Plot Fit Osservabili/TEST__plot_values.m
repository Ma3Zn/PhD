close all
clear
clc

data = "Data/28-2-25";

load(strcat(data, "/exp_values_4_tstep.txt"));
load(strcat(data, "/exp_values_5_tstep.txt"));
load(strcat(data, "/exp_values_6_tstep.txt"));
load(strcat(data, "/folding_factors.txt"));

data_range = 1:5;

exp_values_4_tstep = exp_values_4_tstep(:, data_range);
exp_values_5_tstep = exp_values_5_tstep(:, data_range);
exp_values_6_tstep = exp_values_6_tstep(:, data_range);
folding_factors = folding_factors(data_range);

[n_pow, ~] = size(exp_values_4_tstep);

for i = 1:n_pow
    figure()
    hold on
    grid on

    plot(folding_factors, exp_values_4_tstep(i,:), '-x');
    plot(folding_factors, exp_values_5_tstep(i,:), '-x');
    plot(folding_factors, exp_values_6_tstep(i,:), '-x');

    legend('4_tstep', '5_tstep', '6_tstep');
end
