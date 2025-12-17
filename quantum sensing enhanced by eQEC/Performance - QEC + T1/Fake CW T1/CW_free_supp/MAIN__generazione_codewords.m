clear
clc

%   SCRIPT per la creazione delle codewords per un sistema ad interazioni
%   in competizioni contro pure dephasing
    
%   Dimensione del sistema
    dim = 6;

%   Teniamo traccia della durata temporale della ricerca
    tic
    
%   Generiamo le codewords relative alla dimensione del qudit attuale
    [ew0, ew1] = genera_codewords(dim, 1e-10);
        
%   Aumentiamo il contatore
    j = j + 1;
    
%   Teniamo traccia della durata temporale della ricerca
    toc

%   Salviamo le codewords su file
    scrivi_codewords(ew0, ew1, dim);