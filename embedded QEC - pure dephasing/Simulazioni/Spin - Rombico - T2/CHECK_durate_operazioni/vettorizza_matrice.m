function [v] = vettorizza_matrice(rho)

    n = length(rho);
    v = sparse(zeros(n^2,1));

    for i = 1:n
        i_v = sparse(zeros(n,1));
        i_v(i) = 1;

        for j = 1:n
            j_v = sparse(zeros(n,1));
            j_v(j) = 1;

            v = v + rho(i,j)*kron(i_v, j_v);
        end
    end

    if (~issparse(rho))
        v = full(v);
    end

end