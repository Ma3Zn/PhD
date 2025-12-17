function E = operatori_errori(dim)
%   Funziona che ritorna gli operatori derivanti dalla tomografia del
%   processo di evoluzione idle del sistema

    nome = strcat('./Tomografia/tom_',num2str(dim),'_lv.mat');
    load(nome,'E');
end