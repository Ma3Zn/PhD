clear
clc

sim = genera_op_QEC(6,1);

[gl_c, gl_r] = genera_gate_logico(sim,pi,1);

ang = pi;

cp = diag([1 1 1 exp(-1i*pi)]);

psi_in_sr = [1 1 1 1].'/2;
psi_fin_sr = cp * psi_in_sr;

[sim.rho, psi_in_r] = crea_condizioni_iniziali(sim, [1/sqrt(2) 1/sqrt(2)], [1/sqrt(2) 1/sqrt(2)], 1);

psi_fin_r = pulisci_matrice(sim.CB_ridotta' * (gl_r * gl_r * psi_in_r), 1e3*eps)

%   Da finire ma dovrei aver capito!!