figure(1)
hold on

genera_input_1_2(5e4, pi, [1 0], [1 0], 10);

for dim = 4:12
%%  CU
    sim = genera_op_QEC(dim, 1);

    CU = genera_CU(4);
    log_CU = pulisci_matrice(logm(full(CU)), 1e3*eps);

%     d = diag(kron(sim.CB, eye((2*dim+1)/2)) * log_CU * kron(sim.CB', eye((2*dim+1)/2)));

%     err((2*dim+1)/2) = norm(d);

    [gl_c, gl_r] = genera_gate_logico(sim,pi,1);
    
    log_gl = pulisci_matrice(logm(full(gl_c)), 1e3*eps);

    d = diag(sim.CB_completa * log_gl * sim.CB_completa');

    log_gl = pulisci_matrice(logm(full(gl_r)), 1e3*eps);

    d = diag(sim.CB_ridotta * log_gl * sim.CB_ridotta');


    err_gl((2*dim+1)/2) = norm(d);

    for k = 1:(2*dim+1)/2
        log_R = pulisci_matrice(logm(full(sim.R(:,:,k))), 1e3*eps);

        d = diag(sim.CB * log_R * sim.CB');

        err_R(k) = norm(d);
    end
    plot(err_R);
end

plot(err);

figure(2)
plot(err_gl);