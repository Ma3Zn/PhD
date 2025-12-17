function [op] = genera_superoperatore_sx(A)
%   A puo' essere sia una matrice piena che una matrice sparsa. Se sparsa,
%   il risultato è una matrice sparsa.
    [n, ~] = size(A);

    if (issparse(A))
        op = kron(A, speye(n));
    else
        op = kron(A, eye(n));
    end
end