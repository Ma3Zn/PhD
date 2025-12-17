clear
clc

%   Fissiamo la dimensione del sistema
dim = 6;

%   Path to W
path = strcat("../Grafico Performance + T1/par/", num2str(dim), "/W.txt");
load(path);

%%  Generazione stato iniziale

%   Generiamo uno stato eccitato
v = sparse(dim,1);
idx = 5;
v(idx) = 1;
% v = [0 1 1 0 0 1].'/sqrt(3);
rho = v*v.';

%%  Lindblad con ratei

%   Fattore di scala --> -1/log(rho_rate(idx,idx))
sT = 0.1114;

%   rateo di forza di decoerenza
T = 1;

%   Tempo di evoluzione temporale
t = 1;

%   Costruiamo il superoperatore incoerente
sop_inc_rate = costruisci_sop_inc_rate(W, dim);

%   Scaliamo il superoperatore incoerente
sop_inc_rate = (t/T) * sT * sop_inc_rate;

%   Costruiamo il superoperatore di evoluzione (T1)
sop_evol_rate = fastExpm(sop_inc_rate);

%   Facciamo evolvere il sistema
rho_rate = devettorizza_matrice(sop_evol_rate * vettorizza_matrice(rho))

%%  Calcolo popolazione 0L ed 1L

supp0 = [2 3 6];
supp1 = [1 4 5];

pop_0 = sum(diag(rho_rate(supp0, supp0)))
pop_1 = sum(diag(rho_rate(supp1, supp1)))