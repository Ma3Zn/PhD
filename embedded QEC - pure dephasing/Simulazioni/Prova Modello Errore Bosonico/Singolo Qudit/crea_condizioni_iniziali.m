function [rho, psi_in] = crea_condizioni_iniziali(sim, amp)
%   Funzione che dati i parametri della simulazione e le ampiezze logiche
%   dei due qudit crea le condizioni iniziali del sistema totale (rho,
%   matrice di densità vettorizzata) e del sistema ridotto (psi_in,
%   vettore di stato per il calcolo della fidelity finale)

%   Creiamo i vettori di stato iniziale (in base computazionale) dei due
%   qudit logici
    psi_in_q = amp(1) * sim.CB(:,1) + amp(2) * sim.CB(:,sim.dim_q/2+1);

%   Creiamo lo stato iniziale delle ancille per la ec [tutte sempre in |0>]
    psi_in_a    = sparse(zeros(sim.dim_a, 1));
    psi_in_a(1) = 1;

%   Creaimo il vettore di stato iniziale del sistema ridotto
    psi_in = psi_in_q;

%   Creiamo il vettore di stato del sistema completo
    psi_in_completo = kron(psi_in_q, psi_in_a);

%   Creiamo e vettoriziamo la matrice di densità iniziale del sistema
%   totale
    rho = vettorizza_matrice(psi_in_completo * psi_in_completo');
end