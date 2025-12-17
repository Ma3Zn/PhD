function [rho_par] = partial_trace_1(rho, dim_0, dim_1)

%   vedi: https://physics.stackexchange.com/questions/179671/how-to-take-partial-trace

    rho_par = sparse(zeros(dim_0, dim_0));
    I = speye(dim_0);

    for i = 1:dim_1
        e_i    = sparse(zeros(dim_1,1));
        e_i(i) = 1;
        
        mat = kron(I,e_i);

        rho_par = rho_par + mat' * rho * mat;
    end

    if (~issparse(rho))
        rho_par = full(rho_par);
    end
end