function [overlap, err] = calcola_overlap(sim, psi)
%   Funzione che calcola l'overlap del nostro stato finale con il ground
%   state effettivo dell'Hamiltoniana
    
    J = sim.J;

%   Costruiamo l'Hamiltoniana di Heisenberg all'istante finale
    H = J.z*sim.H_heis.Sz_odd  + J.x*sim.H_heis.Sx_odd  + J.y*sim.H_heis.Sy_odd ...
      + J.z*sim.H_heis.Sz_even + J.x*sim.H_heis.Sx_even + J.y*sim.H_heis.Sy_even;

% % % %   DBG
% % %     H = H + H_heis.i;

%   Calcoliamo il ground state dell'Hamiltoniana
    [V,e] = eig(H);

%   Nel caso in cui il ground state sia degenere recuperiamo tutto
%   l'autospazio relativo
    e_min = min(diag(e));
    idxs  = find(abs(diag(e) - e_min) < 1e-12 );

%   Inizializiamo l'overlap
    overlap = 0;

%   Proiettiamo lo stato nell'autospazio del ground state
    for idx = idxs'
        overlap = overlap + abs(V(:,idx)' * psi)^2;
    end

%   Calcoliamo l'errore relativo in norma due dell'energia dello stato
%   finale
    err = norm(e_min - psi' * H * psi) / norm(e_min);
end