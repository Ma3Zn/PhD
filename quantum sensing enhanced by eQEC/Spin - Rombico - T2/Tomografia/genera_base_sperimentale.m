function A = genera_base_sperimentale(E, dim)
%   Funzione che genera una base di operatori partendo dall'operatore
%   fornito e completandola a base ON tramite procedura di GS con gli
%   operatori |n><m|

    A = cell(1,dim);

%   Condizione iniziale
    A{1} = E;

    idx = 1;

%   Cicliamo sugli altri operatori della base
    for n = 1:dim-1
        for m = 1:dim
%           Creiamo l'operatore attuale
            op = zeros(dim, dim);
            op(n,m) = 1;

%           Ortogonalizziamo tramite GS
            A{idx+1} = matrix_gs(A, op, idx);

%           Aggiorniamo il contatore degli operatori inseriti
            idx = idx + 1;
        end
    end
    
%   Dobbiamo generare ancora dim-1 operatori
    for m = dim:-1:2
%       Creiamo l'operatore attuale
        op = zeros(dim, dim);
        op(dim, m) = 1;

%       Ortogonalizziamo tramite GS
        A{idx+1} = matrix_gs(A, op, idx);

%       Aggiorniamo il contatore degli operatori inseriti
        idx = idx + 1;
    end
end