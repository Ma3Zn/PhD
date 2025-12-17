clear
clc

h_bar = 6.582119569000000e-04;
mu_b  = 0.057883826350000e+00;
g = 2;

%%  Operatori di spin

spin = 3/2;
dim = 2*spin + 1;

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

k = 1
e_vals = linspace(-0.03,-0.01, 1e3);

for e = e_vals
    
    Ham = d * Sz * Sz + e * (Sx * Sx - Sy * Sy) + g * mu_b * B * (Sz * cos(a) + Sx * sin(a));
    Ham = Ham / (spin * (2*spin - 1));
    [V, H0] = eig(Ham);
    
    H0 = H0(end:-1:1, end:-1:1);
    V  =  V(   :    , end:-1:1);
    
    % %   Riordiniamo gli autovettori per essere consistenti con l'ordinamento di
    % %   Sz
    % V = V(:, end:-1:1);
    
    Sxr = V' * Sx * V;
    Syr = V' * Sy * V;
    Szr = V' * Sz * V;
    
    m_Sx(k) = min(min(abs(Sxr(abs(Sxr)>1e-10))));
    m_Sz(k) = min(min(abs(Szr(abs(Szr)>1e-10))));

    k = k+1;
end

semilogy(e_vals, m_Sx, e_vals, m_Sz)
grid on