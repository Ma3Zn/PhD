clear
clc

load tom_3_2_SEQ_CU_MOD_MOD_L.mat

%%

N = 1e3;
err = zeros(N,1);

sz = size(E);

% Test influenza E_1;
% E(:,:,1) = zeros(in.dim);

for i = 1:N
    psi_0 = random_superposition(4);
    psi_1 = random_superposition(2);
    
    psi = kron(psi_0, psi_1);

    rho_init = psi * psi';

    rho_fin_corretta = evoluzione_QO(rho_init, in);

    rho_fin_tom = zeros(8,8);

    for j = 1:sz(3)
        if (norm(E(:,:,j), 'fro') > 1e-2)
            rho_fin_tom = rho_fin_tom + pulisci_matrice(E(:,:,j), 1e-4) * rho_init * pulisci_matrice(E(:,:,j)', 1e-4);
        end
    end

    err(i) = norm(in.U * rho_fin_tom * in.U' - rho_fin_corretta, 'fro');
end

semilogy(err)