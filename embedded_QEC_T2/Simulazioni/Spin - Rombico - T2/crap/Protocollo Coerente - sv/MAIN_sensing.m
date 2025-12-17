close all
clear
clc

%   SCRIPT per la simulazione di una rotazione logica indotta dal campo
%   perpendicolare esterno grazie ad un impulso oscillante lungo Z che
%   rende degenere i livelli opportuni

%   Generazione parametri sistema
sys = genera_sistema();

%   Valori di Bx da testare
Bx_vals = 2e-3;

%   Contatore
i = 1;

%   Cicliamo su vari valori del Bx esterno
for Bx = Bx_vals
    tic

%   Aggiorniamo il valore di Bx da misurarefor
    sys.Bx = Bx;

%   Recuperiamo la frequenza di oscillazione indotta da Bx
    [sim, simfreq(i)] = calcola_frequenza_oscillazione(sys, Bx);

%   Aggiorniamo il contatore
    i = i + 1;

    toc
end

save ws.mat