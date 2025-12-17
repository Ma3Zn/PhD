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

Sxr = V' * Sx * V
Syr = V' * Sy * V
Szr = V' * Sz * V

%%  Scrittura su file di GAMMA
save G.txt G -ascii

save ../Tomografia/par/G.txt G -ascii

save ../CW_T2/'Performance Codewords ATA'/par/Gamma/G.txt G -ascii
save ../CW_T2_free_supp/'Performance Codewords ATA'/par/Gamma/G.txt G -ascii