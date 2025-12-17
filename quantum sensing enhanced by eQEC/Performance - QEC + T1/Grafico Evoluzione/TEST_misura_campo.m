clear
clc

% Script per simulazioni Monte Carlo dell'esecuzione ideale del protocollo
% di sensing con QEC


% Spin del sistema
S = 5/2;

% Valore del campo da misurare [ T ]
Bx = 2e-9;
delta = 1e-10;

% Tempo di coerenza del sistema
T2 = 1e4;
% Generiamo il sistema
tic
sim = genera_parametri_simulazione(S, Bx, T2, 2);
toc

% Eseguiamo le simulazioni del protocollo
tic
[T_es(1), T_mis(1)] = esegui_protocollo_sensing(sim, 0.7);
toc

tic
sim = genera_parametri_simulazione(S, Bx - Bx*delta, T2, 2);
[T_es(2), T_mis(2)] = esegui_protocollo_sensing(sim, 0.7);
toc

% Approssimiamo la derivata dT/dB
dSdB = abs( (T_mis(1) - T_mis(2)) / (Bx*delta) );

%%

% Valutiamo l'efficacia del protocollo
[ Bx, T_es(1)*1e-6, 1/dSdB ]

%%  Salviamo i risultati

save(strcat('5_2_1e4_0.7_2--', string(datetime), '__.mat'))