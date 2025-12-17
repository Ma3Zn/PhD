function U = genera_operatore_rf(imp, t)
%   Funzione per la costruzione dell'operatore di cambio di base tra
%   interaction picture e rotating frame

%   Recuperiamo la lunghezza della sequenza di impulsi da simulare  
    n_imp = length(imp);

%   Allochiamo lo spazio per il logaritmo della matrice da esponenziare
    lU = zeros(4,4);

    for k = 1:n_imp
%       Recuperiamo l'impulso attuale
        p = imp{k};
        
        tmp = zeros(4,4);
        tmp(p.lv(1), p.lv(1)) = 1;
        tmp(p.lv(2), p.lv(2)) = 1;

%       Costruiamo l'operatore di cui calcolare l'esponenziale
        lU = lU + (p.w * t * tmp + p.Bz/p.w * sin(p.w*t) * diag(diag(p.Mz)));
    end

%   Eseguiamo l'esponenziazione
    U = diag(exp(1i * diag(lU)));
end
