function [y] = gs(A,v)
% Funzione che esegue il processo di ortogonalizzazione di Gram-schmidt su
% v dati i vettori colonna di A
    [n,m] = size(A);

    y   = zeros(n,1);
    sum = zeros(n,1);
    u   = zeros(n,1);

    for i = 1:m
        u = A(:,i);

        sum = sum + (v'*u)/(u'*u).*u;
    end

    y = v - sum;
    y = y / sqrt((y'*y));
end