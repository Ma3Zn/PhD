function [B] = pulisci_matrice(A, tol)
%   Secondo me si puo' migliorare ulteriormente ma ora non ne ho il tempo

    [R,C,val] = find(A);

    R_val = real(val);
    I_val = imag(val);

    [n,m] = size(A);

    B      = sparse(n,m);
    B_real = sparse(n,m);
    B_imag = sparse(n,m);

    idxs = find(abs(R_val) - tol > 0);
    for i = idxs'
        B_real(R(i), C(i)) = real(val(i));
    end

    idxs = find(abs(I_val) - tol > 0);
    for i = idxs'
        B_imag(R(i), C(i)) = imag(val(i));
    end

    B = B_real + 1i*B_imag;
end