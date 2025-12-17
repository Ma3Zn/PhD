function [rho] = devettorizza_matrice(v)
    n = length(v);
    n = n^0.5;

    rho = sparse(zeros(n,n));

    for i = 1:n
        i_v = sparse(zeros(n,1));
        i_v(i) = 1;

        for j = 1:n
            j_v = sparse(zeros(n,1));
            j_v(j) = 1;

            rho(i,j) = kron(i_v,j_v)' * v;
        end
    end

    if (~issparse(v))
        rho = full(rho);
    end
end