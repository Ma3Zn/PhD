clear
clc

figure(1)
hold on

genera_input_1_2

dim = 8;
% for dim = 4:2:12
%%  CU
    sim = genera_input_QEC(dim);

    log_CU = pulisci_matrice(logm(full(sim.CU)), 1e3*eps);

    d = diag(kron(sim.CB, eye(dim/2)) * log_CU * kron(sim.CB', eye(dim/2)));

    log_CU = pulisci_matrice(kron(sim.CB, eye(dim/2)) * log_CU * kron(sim.CB', eye(dim/2)), 1e-12);
    [r,c,v] = find(log_CU*i);
    theta = abs(v);

% % %     full(sparse(r,c,theta))
% % %     full(sparse(r,c,angle(v./theta)))

    err(dim/2) = norm(d);

%% gate logico
    [gl, glr] = genera_gate_logico(sim, ang);
    sim.gate_logico = full(gl);
    sim.gate_logico_ridotto = full(glr);

    log_gl = pulisci_matrice(logm(full(sim.gate_logico)), 1e3*eps);

    d = diag(log_gl);

    log_glr = full(pulisci_matrice(logm(sim.gate_logico_ridotto), 1e-12))
    [r,c,v] = find(log_glr*i);
    theta = abs(v);

    full(sparse(r,c,theta))
    full(sparse(r,c,angle(v./theta)))
    
    err_gl(dim/2) = norm(d);

    for k = 2:dim/2
        log_R = pulisci_matrice(logm(full(sim.R{k})), 1e3*eps);

        d = diag(log_R);

        [r,c,v] = find(log_R*i);
        theta = abs(v);

% % %         full(sparse(r,c,theta))
% % %         full(sparse(r,c,angle(v./theta)))

        err_R(k) = norm(d);
    end
    plot(err_R);
% end

plot(err);

figure(2)
plot(err_gl);