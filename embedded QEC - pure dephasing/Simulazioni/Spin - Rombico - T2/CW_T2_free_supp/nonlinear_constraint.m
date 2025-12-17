function [cineq,ceq] = nonlinear_constraint(x)
%   Funzione che calcola il valore dei vincoli non lineari dato il vettore
%   x in input, tale è un vettore di dimensione 2*d contenente gli elementi
%   della codifica di |0_L> nelle prime d entrate e quelli della codifica
%   di |1_L> nelle restanti d.
%   Necessità di una struttura globale ERR contenete lo spacing delle
%   codewords.
%   Abbiamo solo vincoli di uguaglianza e dunque cineq = []

    global ERR;

%   Recuperiamo lo spacing della codifica delle codewords
    s = ERR.spacing;

%   Recuperiamo la dimensione del sistema
    dim = ERR.dim;

%   Settiamo vuoti i vincoli di inequaglianza
    cineq = [];

%   Inizializziamo i vincoli di uguaglianza
    ceq = [];

%   Calcoliamo i vincoli di normalizzazione
    ceq(1) = x(1:dim)' * x(1:dim) - 1;
    ceq(2) = x(dim+1:2*dim)' * x(dim+1:2*dim) - 1;

%   Teniamo traccia del numero di vincoli di uguaglianza
    n_v = 2;

%   Calcoliamo i vincoli di spacing 
    for i = 0:dim-1

%       Calcoliamo i bound su j
        lb = dim + max(0,i-(s-1));
        ub = dim + min(dim-1,i+(s-1));

        for j = lb:ub
%           Calcoliamo il vincolo di spacing
            ceq(n_v + 1) = x(i+1) * x(j+1);

%           Penality
            ceq(n_v + 1) = 1e1 * ceq(n_v + 1);

%           Aggiorniamo il numero di vincoli
            n_v = n_v + 1;
        end
    end

%   Includiamo le condizioni di Knill-Laflamme come vincoli e dunque
%   trasformiamo il problema da un problema di ottimizzazione vincolato ad
%   un (forse) più semplice problema di soddisfacimento di vincoli
    KL = KL_I(x);

%   Recuperiamo il numero di condizioni di Knill-Laflamme
    n_kl = length(KL);

%   Aggiungiamo le condizioni di Knill-Laflamme al vettore dei vincoli
    ceq(n_v+1:n_v + n_kl) = KL;

%   Aggiorniamo il numero di vincoli di uguaglianza
    n_v = n_v + n_kl;

% % % %   Aggiungiamo la forzatura dei supporti delle codewords
% % %     for j = ERR.idx0
% % %         ceq(n_v + 1) = x(j);
% % %         n_v = n_v + 1;
% % %     end
% % % 
% % %     for j = ERR.idx1
% % %         ceq(n_v + 1) = x(j);
% % %         n_v = n_v + 1;
% % %     end
end