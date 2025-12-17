clear
clc

%   Fissiamo la dimensione del sistema
dim = 6;

%   Leggiamo i ratei di T1
G = load("G.txt");
  
% fattore di scala consistente col T2 del caso FM
% G = G * 0.7;
% così è forse un po' troppo veloce il T2, è comunque l'ultima cosa da
% sistemare
G = G;

%   Riduciamo le dimensioni di G
G = G(1:dim, 1:dim);
g = G;

%   Ricaviamo gamma da GAMMA
for mu=1:dim
    for nu=1:dim
        g(mu,nu) = 2*G(mu,nu)-G(mu,mu)-G(nu,nu);
    end
end

%   Generiamo uno stato eccitato
v = sparse(dim,1);
v(1) = 1/sqrt(2);
v(2) = 1/sqrt(2);
rho  = v*v.';

%   rateo di forza di T2
T2 = 1e5;

%   Tempo di evoluzione temporale
t = 2e5;

%   Costruiamo il superoperatore incoerente
sop_inc_T2 = costruisci_sop_inc_T2(g, dim);

%   Scaliamo il superoperatore incoerente
sop_inc_T2 = (t/T2) * sop_inc_T2;

%   Costruiamo il superoperatore di evoluzione (T1)
sop_evol_T2 = fastExpm(sop_inc_T2);

%   Facciamo evolvere il sistema
rho_T2 = devettorizza_matrice(sop_evol_T2 * vettorizza_matrice(rho))

