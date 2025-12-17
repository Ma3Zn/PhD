function [cineq,ceq] = nonlinear_constraint(x)
%   Funzione che calcola il valore dei vincoli non lineari dato il vettore
%   x in input, tale è un vettore di dimensione 2*d contenente gli elementi
%   della codifica di |0_L> nelle prime d entrate e quelli della codifica
%   di |1_L> nelle restanti d.
%   Necessità di una struttura globale ERR contenete lo spacing delle
%   codewords.
%   Abbiamo solo vincoli di uguaglianza e dunque cineq = []

    global ERR;

%   Recuperiamo la dimensione del sistema
%   Poichè forziamo i supporti e dunque x contiene solamente gli lementi
%   non nulli, nel nostro caso la dimensione corretta è dim/2
    dim = ERR.dim/2;

%   Settiamo vuoti i vincoli di inequaglianza
    cineq = [];

%   Inizializziamo i vincoli di uguaglianza
    ceq = [];

%   Calcoliamo i vincoli di normalizzazione
    ceq(1) = x(1:dim)' * x(1:dim) - 1;
    ceq(2) = x(dim+1:2*dim)' * x(dim+1:2*dim) - 1;

%   Teniamo traccia del numero di vincoli di uguaglianza
    n_v = 2;

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

% ATTENZIONE::  Dobbiamo rimovere tutta la parte immaginaria dai vincoli in
%               quanto prima fa ottimizzazione reale

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