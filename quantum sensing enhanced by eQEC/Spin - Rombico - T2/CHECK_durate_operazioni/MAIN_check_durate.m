clear
clc

h_bar = 6.582119569000000e-04;
mu_b  = 0.057883826350000e+00;
g = 2;

%%  Parametri sistema

spin = 7/2;
dim = 2*spin + 1;

sim = genera_input_QEC(dim);

%%  Operatori di spin
Sx = sop(spin, 'x');
Sy = sop(spin, 'y');
Sz = sop(spin, 'z');

%%  Parametri del sistema

d = -0.1;
e_max = d / 3;
e = -0.03;

% PLOT autovalori
% d = 1/3 * e; % S = 3/2
% d = 0.23 * e; % S = 7/2

B = 0.35; % T
a = 0;

Ham = d * Sz * Sz + e * (Sx * Sx - Sy * Sy) + g * mu_b * B * (Sz * cos(a) + Sx * sin(a));
Ham = Ham / (spin * (2*spin - 1));
[V, H0] = eig(Ham);

H0 = H0(end:-1:1, end:-1:1);
V  =  V(   :    , end:-1:1);

% %   Riordiniamo gli autovettori per essere consistenti con l'ordinamento di
% %   Sz
% V = V(:, end:-1:1);

%%  Scrittura di GAMMA

Sxr = diag(V' * Sx * V);
Szr = diag(V' * Sz * V);

G = zeros(dim, dim);

for n = 1:dim
    for m = 1:dim
        G(m,n) = Sxr(m)*Sxr(n) + Sxr(m)*Szr(n) + Szr(m)*Sxr(n) + Szr(m)*Szr(n);
    end
end

Sxr = V' * Sx * V;
Syr = V' * Sy * V;
Szr = V' * Sz * V;

min(min(abs(Sxr(abs(Sxr)>1e-10))))
min(min(abs(Szr(abs(Szr)>1e-10))))

%%  Costruzione Generatore Rotazione Logica
% Scriviamo il gate logico in base computazionale
[~, gl_ridotto] = genera_gate_logico(sim, [pi/2, 0]);

% Calcoliamo il generatore del gate logico 
log_G  = -1i * logm(full(gl_ridotto));

%%  Calcolo campo opportuno
t_fix = 90; % ns

B = h_bar / (g * mu_b * t_fix) ...
  * max(max( abs(log_G(abs(log_G) > 1e-12)) ./ abs(Sxr(abs(log_G) > 1e-12)) ))

%%  Calcolo tempo opportuno

B_fix = 0.005; % mT --> 50G

T = h_bar / (g * mu_b * B_fix) ...
  * max(max( abs(log_G(abs(log_G) > 1e-12)) ./ abs(Sxr(abs(log_G) > 1e-12)) ))