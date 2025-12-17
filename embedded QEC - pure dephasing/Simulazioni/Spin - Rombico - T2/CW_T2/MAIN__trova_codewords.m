clear
clc

%   SCRIPT per la creazione delle codewords per un sistema ad interazioni
%   in competizioni

%   Supporti delle codifiche
global supp0;
global supp1;

% % % %   Definiamo i suporti della codifica (4 lv)
% % % supp0 = [2 3];
% % % supp1 = [1 4];

% % % %   Definiamo i suporti della codifica (6 lv)
% % % supp0 = [2 3 6];
% % % supp1 = [1 4 5];

%  Definiamo i suporti della codifica (8 lv)
supp0 = [2 3 6 7];
supp1 = [1 4 5 8];

% % % %   Definiamo i suporti della codifica (10 lv)
% % % supp0 = [2 3 6 7 10];
% % % supp1 = [1 4 5 8  9];

%   Dimensione del sistema
dim = 8;
    
%   Generiamo le codewords relative alla dimensione del qudit attuale
[ew0, ew1, pass] = genera_codewords(dim, 1e-14);

%   Salviamo le codewords su file
scrivi_codewords(ew0, ew1);