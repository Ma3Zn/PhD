% Script per la process tomography del CU nella procedura di QEC. Gli
% operatori della quantum operation ottenuta vengono poi decompositi sulla
% base degli operatori di Steevens

clear
clc

% dimensione del Qudit
dim = 4;

%   Generiamo un parpool per l'esecuzione in parallelo della tomografia
% parpool(min(dim,8));
fprintf("Terminata creazione parallel pool.\n\nInzio tomografia...\n");

% Teniamo traccia del tempo di esecuzione del programma
tic

% Generiamo gli input per il nostro algoritmo
in = genera_input_PT(dim);

% Creiamo l'operatore di evoluzione che caratterizza la nostra quantum
% operation
in.evol = genera_operatore_evoluzione_QO(in);

% Eseguiamo la process tomography della nostra quantum operation
[E, CHI] = process_tomography(in);

% Riordiniamo gli operatori in base alla loro norma
E = riordina_operatori(E, dim^2);

nome = strcat("tom_", num2str(dim), "_lv.mat");
save(nome);

nome_ = strcat("../CW_T2/Tomografia/",nome);
save(nome_);

nome_ = strcat("../CW_T2_free_supp/Tomografia/",nome);
save(nome_);

toc