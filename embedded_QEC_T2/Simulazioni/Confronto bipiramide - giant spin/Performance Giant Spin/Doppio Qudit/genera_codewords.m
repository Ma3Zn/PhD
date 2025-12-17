function M = genera_codewords(dim)
%   Funzione che datala dimensione dello spazio e il tempo a cui sono state
%   ottimizzateritorna le opportune codewords e gli indii dei loro supporti

    radice = strcat('./par/', num2str(dim), '/');

%   Apriamo gli stream su file

    nome_0 = strcat(radice, 'cw_0.bin');
    nome_1 = strcat(radice, 'cw_1.bin');

    fileID_0 = fopen(nome_0, 'r');
    fileID_1 = fopen(nome_1, 'r');

%   Leggiamo le due codifiche
    ew0 = fread(fileID_0, 'double');
    ew1 = fread(fileID_1, 'double');

%   Formattiamo correttamente le due codifiche
    M  = [reshape(ew0, dim, dim/2), reshape(ew1, dim, dim/2)];

%   Chiudiamo gli stream su file
    fclose(fileID_0);
    fclose(fileID_1);
end