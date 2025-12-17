function [op] = genera_superoperatore_dx(A)
%   A puo' essere sia una matrice piena che una matrice sparsa. Se sparsa,
%   il risultato è una matrice sparsa.
    [n, ~] = size(A);

    if (issparse(A))
        op = kron(speye(n),A.');
    else
        op = kron(eye(n),A.');
    end
end