spin = 1/2

Sz = diag(-1/2:1/2);

E_0 = Sz;

global alpha
global beta

alpha = 1/sqrt(2);
beta  = 1i/sqrt(2);

psi_in = alpha * [1 0]' + beta * [0 1]';
rho_in = psi_in*psi_in';

rho = vettorizza_matrice(rho_in);

%% Genriamo l'operatore logico da eseguire

ang(1) = 6.842619328173;
ang(2) = 4.917465198531;

theta = rand() * 2*pi;
THETA = rand() * 2*pi;

% op = [0,1;-1,0];
% op = [1,0;0,exp(1i*theta)];

% Non è un Rz ma Rz è il nome che nel coice per la simulaizone
% op = [cos(theta),sin(theta)*exp(-1i*THETA);-sin(theta)*exp(1i*THETA),cos(theta)];  

% Per simulare una fase IDLE
op = eye(2);

%% Generiamo le condizioni finali

psi_fin = op * psi_in;

%% Propaghiamo il gate fisico eseguito

global gate_fisico
gate_fisico = op;