function [Y] = matrix_gs(A, E, idx)
% Funzione che esegue il processo di ortogonalizzazione di Gram-schmidt su
% v dati i vettori colonna di A
    [dim, ~] = size(E);

    sum = zeros(dim,dim);

    for i = 1:idx
        u = A{i};

        sum = sum + trace(E'*u)/trace(u'*u) * u;
    end

    Y = E - sum;
    Y = Y / sqrt(trace(Y'*Y));
end