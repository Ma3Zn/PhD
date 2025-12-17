clear
% clc

%   Script per la generazione di una tomografia fasulla: da utilizzare per
%   forzare gli operatori per cui trovare codewords

%   Spin del sistema considerato
spin = 5/2;

%   Dimensione del sistema
dim = 2*spin + 1;

%   Operatori d'errore del sistema
E = cell(1, dim^2);

%   Generiamo l'operatore d'errore desiderato
E{1} = genera_O21(spin);

%   Inizializzazione dei restatnti operatori d'errore
for i = 2:dim^2
    E{i} = zeros(dim);
end

%   Numero di operatori d'errore da considerare
n = 1;

save ../CW_free_supp/Tomografia/tom_6lv.mat