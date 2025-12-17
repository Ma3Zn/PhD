% 1/2 --> Bx = 5.7223e-6
% 3/2 --> Bx = 2.4333e-7
% 5/2 --> Bx = 4.1423e-8

fig = openfig("./images/Bx(dT).fig");

plot([1e-5 1e-3], [5.7223e-6 5.7223e-6], '--', 'LineWidth', 3);
plot([1e-3 1e-1], [2.4333e-7 2.4333e-7], '--', 'LineWidth', 3);
plot([1e-2 1e-0], [4.1423e-8 4.1423e-8], '--', 'LineWidth', 3);