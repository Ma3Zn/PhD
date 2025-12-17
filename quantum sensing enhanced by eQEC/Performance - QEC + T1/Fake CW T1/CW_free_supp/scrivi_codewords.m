function [] = scrivi_codewords(ew0, ew1, dim)
%   Funzione per la scrittura su file delle codewords del sistema

%   Creiamo la cartella realtiva alle codewords ottenute al tempo di
%   ottimizzazione attuale
    dir = strcat('./Codewords/', num2str(dim));
    mkdir(dir);

%   Creiamo il nome del file
    nome0 = strcat(dir,'/cw_0.bin');
    nome1 = strcat(dir,'/cw_1.bin');

%-------------------------------------------------------------------------%
%   Aprima lo stream di output
    fileID_0 = fopen(nome0,'w');

%   Scriviamo le codewords
    fwrite(fileID_0,ew0,'double');

%   Chiudiamo lo stream di output
    fclose(fileID_0);
%-------------------------------------------------------------------------%
%   Aprima lo stream di output
    fileID_1 = fopen(nome1,'w');

%   Scriviamo le codewords
    fwrite(fileID_1,ew1,'double');

%   Chiudiamo lo stream di output
    fclose(fileID_1);
%-------------------------------------------------------------------------%

end