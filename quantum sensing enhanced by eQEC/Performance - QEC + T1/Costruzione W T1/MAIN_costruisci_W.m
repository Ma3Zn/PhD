clear
clc

h_bar = 6.582119569000000e-04;
mu_b  = 0.057883826350000e+00;
g = 2;

%%  Operatori di spin

spin = 5/2;
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

%% Scrittura di W

%  Costruiamo l'operatore alla base del T1 
op = Sx*Sz + Sz*Sx; % O21  --> sf = 0.111352570916444
op = Sx*Sx - Sy*Sy; % O22  --> sf = 0.269102042436527 --> non ci piace perchè ha elementi di Sxr nulli
op = Sy*Sz + Sz*Sy; % O2-1 --> sf = 9.756468937172134e-02
op = (Sx*Sz + Sz*Sx + Sx*Sx - Sy*Sy)/2; % O21 + O22 --> sf = 8.908205673315521e-02

%   Ruotiamo sulla base degli autovettori
op = V' * op * V;

%   Rimuoviamone l'eventuale diagonale
op = op - diag(diag(op));

W = triu(abs(op).^2)';

%% Salviamo W (NON sovrasciviamo G per evitare di fare dei casini)

%   Utilizzo il formato ASCII solo perchè anche G è stata salvata in tale
%   formato. Per una maggiore precisione in futuro potrebbe convenire
%   salvare tutto quanto in formato nativo MATLAB

path = strcat("../Grafico Performance + T1/par/",num2str(dim),"/W.txt");
save(path, "W", "-ascii");