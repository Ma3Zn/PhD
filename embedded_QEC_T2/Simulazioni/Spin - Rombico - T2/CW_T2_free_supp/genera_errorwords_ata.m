function [ew0, ew1] = genera_errorwords_ata(cw)
%   Funzione che genera le errorwords del sistema dato

%   Variabile globale contenente i dati sull'errore del sistema
    global ERR;

%   Inizializziamo le variabili
    ew0 = zeros(ERR.dim, ERR.n);
    ew1 = zeros(ERR.dim, ERR.n);

%   La prima errorword coincide con la codewords
%     ew0 = cw(:,1);
%     ew1 = cw(:,2);

%   Cicliamo sugli errori del sistema
    for j = 1:ERR.n
%-------------------------------------------------------------------------%
        tmp = ERR.E{j} * cw(:,1);
        tmp = tmp / sqrt(tmp'*tmp);

        ew0(:,j) = gs(ew0(:,1 : j - 1), tmp);
%-------------------------------------------------------------------------%
        tmp = ERR.E{j} * cw(:,2);
        tmp = tmp / sqrt(tmp'*tmp);

        ew1(:,j) = gs(ew1(:,1 : j - 1), tmp);
%-------------------------------------------------------------------------%
    end
end