function [] = scrivi_codewords_ata_T2(ew0, ew1, pass, t_opt)
%   Funzione per la scrittura su file delle codewords del sistema

%   Creiamo la cartella realtiva alle codewords ottenute al tempo di
%   ottimizzazione attuale
    dir = './Codewords/';
    mkdir(dir);

%     for j = 1:length(ew0)<<<
%         if pass(j) == true
%           Recuperiamo la dimensione del sistema
            dim = length(ew0{1});

%           Creiamo il nome del file
            nome0 = strcat(dir,'/cw_',num2str(dim),'lv_0.bin');
            nome1 = strcat(dir,'/cw_',num2str(dim),'lv_1.bin');

%-------------------------------------------------------------------------%
%           Aprima lo stream di output
            fileID_0 = fopen(nome0,'w');

%           Scriviamo le codewords
            fwrite(fileID_0,ew0{1},'double');

%           Chiudiamo lo stream di output
            fclose(fileID_0);
%-------------------------------------------------------------------------%
%           Aprima lo stream di output
            fileID_1 = fopen(nome1,'w');

%           Scriviamo le codewords
            fwrite(fileID_1,ew1{1},'double');

%           Chiudiamo lo stream di output
            fclose(fileID_1);
%-------------------------------------------------------------------------%

%         end
%     end
end