close all
clear
clc

data = "Data/2025-03-20";

backend_noise = "/Brisbane";
% backend_noise = "/Fez";
% backend_noise = "/Future";

dim_SbE = "/H_pow_1";

dir = strcat(data, backend_noise, dim_SbE);

load(strcat(dir, "/exp_values_3_tstep.txt"));
load(strcat(dir, "/exp_values_4_tstep.txt"));
load(strcat(dir, "/exp_values_5_tstep.txt"));
load(strcat(dir, "/exp_values_6_tstep.txt"));
load(strcat(dir, "/exp_values_7_tstep.txt"));
load(strcat(dir, "/folding_factors.txt"));

data_range = 1:17;

exp_values_3_tstep = exp_values_3_tstep(:, data_range);
exp_values_4_tstep = exp_values_4_tstep(:, data_range);
exp_values_5_tstep = exp_values_5_tstep(:, data_range);
exp_values_6_tstep = exp_values_6_tstep(:, data_range);
exp_values_7_tstep = exp_values_7_tstep(:, data_range);
folding_factors = folding_factors(data_range);

[n_pow, ~] = size(exp_values_4_tstep);

for i = 1:n_pow
    figure()
    hold on
    grid on


    plot(folding_factors, exp_values_3_tstep(i,:), '-x');
    plot(folding_factors, exp_values_4_tstep(i,:), '-x');
    plot(folding_factors, exp_values_5_tstep(i,:), '-x');
    plot(folding_factors, exp_values_6_tstep(i,:), '-x');
    plot(folding_factors, exp_values_7_tstep(i,:), '-x');

    legend('3_tstep', '4_tstep', '5_tstep', '6_tstep', '7_tstep');
end
