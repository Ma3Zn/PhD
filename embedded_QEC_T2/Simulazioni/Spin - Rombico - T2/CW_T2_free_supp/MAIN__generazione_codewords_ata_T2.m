clear
clc

%   SCRIPT per la creazione delle codewords per un sistema ad interazioni
%   in competizioni contro pure dephasing

%   Settiamo il tempo a cui eseguire l'ottimizzazione
t_opt = [5e-2 5e-6 5e-10 5e-14];
t_opt = 5e-2;

for t = t_opt
%   Output per tenere traccia dell'avanzamento della simulazione
    t

%   Contatore
    j = 1;
    
%   Inizializziamo le variabili che conterrà le codewords+errorwords
    ew0 = cell(1,7);
    ew1 = cell(1,7);
    
%   Inizializziamo la variabile che conterrà la bontà delle codewords
    pass = zeros(1,7);
    
%   Cicliamo sulla dimensione
    for dim = 8
%       Teniamo traccia della durata temporale della ricerca per ogni istanza
        tic

%       Output per tenere traccia dell'avanzamento della simulazione
        dim
    
%       Generiamo le codewords relative alla dimensione del qudit attuale
        [ew0{j}, ew1{j}, pass(j)] = genera_codewords_ata_T2(dim, t, 1e-14);
        
%       Aumentiamo il contatore
        j = j + 1;
    
%       Teniamo traccia della durata temporale della ricerca per ogni istanza
        toc
    end

%   Salviamo le codewords su file
    scrivi_codewords_ata_T2(ew0, ew1, pass, t);
end