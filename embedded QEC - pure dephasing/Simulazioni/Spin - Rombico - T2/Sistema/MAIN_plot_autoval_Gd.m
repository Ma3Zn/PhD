clear
clc

h_bar = 6.582119569000000e-04;
mu_b  = 0.057883826350000e+00;
g = 2;

%%  Operatori di spin

spin = 7/2;
dim = 2*spin + 1;

Sx = sop(spin, 'x');
Sy = sop(spin, 'y');
Sz = sop(spin, 'z');

%%  Parametri del sistema

% meV --> cm^-1 --> GHz
% *8.0657 --> *30

d = 1281 / 1000 / 30 / 8.0657;
e = 294 / 1000 / 30 / 8.0657;

B = 0.35; % T
a = 0;

Ham = d * Sz * Sz + e * (Sx * Sx - Sy * Sy) + g * mu_b * B * (Sz * cos(a) + Sx * sin(a));
[V, H0] = eig(Ham);

H0 = H0(end:-1:1, end:-1:1);
V  =  V(   :    , end:-1:1);

% %   Riordiniamo gli autovettori per essere consistenti con l'ordinamento di
% %   Sz
% V = V(:, end:-1:1);

for i = 1:dim
    for j = 1:dim
        w(i,j) = (H0(i,i) - H0(j,j)) * 8.0657;
    end
end
w

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

max(w);

[val, I] = min(diff(diff>0))

%% PLOT autovalori
aut = diag(H0) * 8.0657;
figure()
hold on
grid on

for i = 1:dim
    plot(linspace(0,1,2),[aut(i), aut(i)], 'k');
end
