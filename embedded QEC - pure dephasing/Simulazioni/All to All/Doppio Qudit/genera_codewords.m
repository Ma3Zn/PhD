function [M, idxs] = genera_codewords(dim)
%   Funzione che datala dimensione dello spazio ritorna le opportune
%   codewords e gli indii dei loro supporti

    idxs = zeros(2, dim/2);

    if dim == 4
        fileID_0 = fopen('./par/data/errw_0_4_lv.bin', 'r');
        fileID_1 = fopen('./par/data/errw_1_4_lv.bin', 'r');

%       Supporto k = 0
        idxs(1,:) = [1 4];
%       Supporto k = 1
        idxs(2,:) = [2 3];
    elseif dim == 6
        fileID_0 = fopen('./par/data/errw_0_6_lv.bin', 'r');
        fileID_1 = fopen('./par/data/errw_1_6_lv.bin', 'r');

%       Supporto k = 0
        idxs(1,:) = [1 3 6];
%       Supporto k = 1
        idxs(2,:) = [2 4 5];
    elseif dim == 8
        fileID_0 = fopen('./par/data/errw_0_8_lv.bin', 'r');
        fileID_1 = fopen('./par/data/errw_1_8_lv.bin', 'r');

%       Supporto k = 0
        idxs(1,:) = [1 2 3 8];
%       Supporto k = 1
        idxs(2,:) = [4 5 6 7];
    elseif dim == 10
        fileID_0 = fopen('./par/data/errw_0_10_lv.bin', 'r');
        fileID_1 = fopen('./par/data/errw_1_10_lv.bin', 'r');

%       Supporto k = 0
        idxs(1,:) = [2 3 5 7 10];
%       Supporto k = 1
        idxs(2,:) = [1 4 6 8 9];
    else
        fileID_0 = fopen('./par/data/errw_0_12_lv.bin', 'r');
        fileID_1 = fopen('./par/data/errw_1_12_lv.bin', 'r');

%       Supporto k = 0
        idxs(1,:) = [2 4 5 7 9 12];
%       Supporto k = 1
        idxs(2,:) = [1 3 6 8 10 11];
    end

    M0 = fread(fileID_0, 'double');
    M1 = fread(fileID_1, 'double');

    M  = [reshape(M0, dim, dim/2), reshape(M1, dim, dim/2)];

    fclose(fileID_0);
    fclose(fileID_1);
end