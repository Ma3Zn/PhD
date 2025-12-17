function M = genera_codewords(dim)
%   Funzione che datala dimensione dello spazio e il tempo a cui sono state
%   ottimizzateritorna le opportune codewords e gli indii dei loro supporti

    nome_0 = strcat('../Codewords/cw_',num2str(dim),'lv_0.bin');
    nome_1 = strcat('../Codewords/cw_',num2str(dim),'lv_1.bin');

    fileID_0 = fopen(nome_0, 'r');
    fileID_1 = fopen(nome_1, 'r');

    ew0 = fread(fileID_0, 'double');
    ew1 = fread(fileID_1, 'double');

    M  = [reshape(ew0,dim,dim/2), reshape(ew1,dim,dim/2)];

    fclose(fileID_0);
    fclose(fileID_1);
end