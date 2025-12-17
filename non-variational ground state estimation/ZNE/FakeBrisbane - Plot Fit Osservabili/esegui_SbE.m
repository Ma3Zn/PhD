function e_SbE = esegui_SbE(H, S)
%   Funzione che risolve il problema generalizzato agli autovalori relativo
%   alla SbE avente per matrici H ed S. Ritorna il più piccolo autovalore
%   trovato

%   Risolviamo il problema agli autovalori generalizzato Hv = eSv con
%   le routine di MATLAB, dettagli su come risolvere tale problema
%   possono essere trovati in un qualsiasi testo di analisi numerica o
%   anche nell'articolo "Quantum Power Method by a Superposition of
%   Time-Evolved States"
    [~, D] = eig(H, S);

%   Recuperiamo l'autovettore relativo al minimo autovalore del
%   problema generalizzato
    e_SbE = min(diag(D));
end