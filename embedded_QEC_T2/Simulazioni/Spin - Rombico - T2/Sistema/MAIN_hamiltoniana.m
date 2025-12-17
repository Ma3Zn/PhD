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
e = -0.01;

% PLOT autovalori
% d = 1/3 * e; % S = 3/2
% d = 0.23 * e; % S = 7/2

B = 1.05; % T
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

min(min(abs(Sxr(abs(Sxr)>1e-10))));
min(min(abs(Szr(abs(Szr)>1e-10))));

%% GAP
k = 1;
for i = 1:dim
    for j = i+1:dim
        w(k) = (H0(i,i) - H0(j,j)) * 242.9;
        k = k+1;
    end
end

%% Diff gap
k = 1;
for i = 1:length(w)
    for j = 1:length(w)
        if (i ~= j)
            diff(k) = abs(w(i) - w(j));
            k = k+1;
        end
    end
end

min(w)

% [val, I] = min(diff(diff>0));

%%  Scrittura su file di GAMMA
save G.txt G -ascii

save ../Tomografia/par/G.txt G -ascii
save ../CW_T2/par/G.txt G -ascii
save ../CW_T2/'Performance Codewords ATA'/par/Gamma/G.txt G -ascii
save ../CW_T2_free_supp/'Performance Codewords ATA'/par/Gamma/G.txt G -ascii