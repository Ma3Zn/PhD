function exp_val_ZNE = calcola_valori_aspettazione_ZNE(exp_val, ff)
%   Funzione che passata una matrice di valori di aspettazione, per ogni
%   riga della matrice ne calcola la ZNE tramite un modello a somma di due
%   esponenziali

%   Recuperiamo le dimensioni della matrice
    [n_pow, ~] = size(exp_val);

    for i = 1:n_pow
%       Eseguiamo il fit dei dati
        fitted_model = fit(ff', exp_val(i,:)', 'exp2');

% % % %       Eseguiamo l'estrapolazione a rumore zero -- singolo esponenziale
% % %         exp_val_ZNE(i) = fitted_model.a;

%       Eseguiamo l'estrapolazione a rumore zero -- doppio esponenziale
        exp_val_ZNE(i) = fitted_model.a + fitted_model.c;

% % % %       Eseguiamo l'estrapolazione a rumore zero -- LINEARE
% % %         exp_val_ZNE(i) = fitted_model.p2;

% % % %       Eseguiamo l'estrapolazione a rumore zero -- quadratica
% % %         exp_val_ZNE(i) = fitted_model.p3;

    end
end
