function [ew0, ew1] = genera_errorwords(cw)
%   Funzione che genera le errorwords del sistema dato

%   Variabile globale contenente i dati sull'errore del sistema
    global ERR;

    global supp0;
    global supp1;

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

        ew0(supp0,j) = gs(ew0(supp0,1 : j - 1), tmp(supp0));
%-------------------------------------------------------------------------%
        tmp = ERR.E{j} * cw(:,2);
        tmp = tmp / sqrt(tmp'*tmp);

        ew1(supp1,j) = gs(ew1(supp1,1 : j - 1), tmp(supp1));
%-------------------------------------------------------------------------%
    end
end