function e_gs = calcola_energia_gs(sim)
%   Funzione per il calcolo dell'energia del ground state dell'Hamiltoniana
%   considerata

%   Funzione che calcola l'overlap del nostro stato finale con il ground
%   state effettivo dell'Hamiltoniana
    J = sim.J;

%   Costruiamo l'Hamiltoniana di Heisenberg all'istante finale
    H = J.z*sim.H_heis.Sz_odd  + J.x*sim.H_heis.Sx_odd  + J.y*sim.H_heis.Sy_odd ...
      + J.z*sim.H_heis.Sz_even + J.x*sim.H_heis.Sx_even + J.y*sim.H_heis.Sy_even;

%   Calcoliamo il ground state dell'Hamiltoniana
    [V,e] = eig(H);

%   Recuperiamo l'energia relativa al ground state
    e_gs = min(diag(e));
end