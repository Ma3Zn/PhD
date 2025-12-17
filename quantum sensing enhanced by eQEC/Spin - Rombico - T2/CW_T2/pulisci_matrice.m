function [B] = pulisci_matrice(A, tol)
    [n,m] = size(A);
    B = zeros(n,m);

    for i = 1:n
        for j = 1:m
            if (abs(real(A(i,j))) > tol)
                B(i,j) = B(i,j) + real(A(i,j));
            end

            if (abs(imag(A(i,j))) > tol)
                B(i,j)= B(i,j) + 1i * imag(A(i,j));
            end
        end
    end
end