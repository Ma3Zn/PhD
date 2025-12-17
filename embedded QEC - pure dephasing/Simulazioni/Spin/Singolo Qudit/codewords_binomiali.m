function [e_00, e_01] = codewords_binomiali(spin)
    dim = 2*spin + 1;

    e_00 = zeros(dim,1);
    e_01 = zeros(dim,1);

    for k = 1:2:dim-1
        j = k+1;

        e_00(j) = sqrt(nchoosek(dim-1,k));
    end

    for k = 0:2:dim-1
        j = k+1;

        e_01(j) = sqrt(nchoosek(dim-1,k));
    end

    e_00 = e_00 / sqrt(2^(dim-2));
    e_01 = e_01 / sqrt(2^(dim-2));
end