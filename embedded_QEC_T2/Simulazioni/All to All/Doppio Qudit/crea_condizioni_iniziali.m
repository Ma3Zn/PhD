function [rho, psi_in] = crea_condizioni_iniziali(sim, amp_c, amp_t, q)
%   Funzione che dati i parametri della simulazione e le ampiezze logiche
%   dei due qudit crea le condizioni iniziali del sistema totale (rho,
%   matrice di densità vettorizzata) e del sistema ridotto (psi_in,
%   vettore di stato per il calcolo della fidelity finale)

%   Creiamo i vettori di stato iniziale (in base computazionale) dei due
%   qudit logici
    psi_in_c = amp_c(1) * sim.CB(:,1) + amp_c(2) * sim.CB(:,sim.dim_q/2+1);
    psi_in_t = amp_c(1) * sim.CB(:,1) + amp_c(2) * sim.CB(:,sim.dim_q/2+1);

%   Creiamo il vettore di stato iniziale (in base computazionale) dello
%   switch logico [sempre in |0L>]
    psi_in_s = sim.CB(:,1);

%   Creiamo lo stato iniziale delle ancille per la ec [tutte sempre in |0>]
    psi_in_a    = sparse(zeros(sim.dim_a, 1));
    psi_in_a(1) = 1;

%   Creaimo il vettore di stato iniziale del sistema ridotto
    psi_in = kron(kron(psi_in_c, psi_in_t), psi_in_s);

%   Creiamo il vettore di stato del sistema completo
    switch q
        case 1
            psi_in_completo = kron(psi_in_c, psi_in_a);
            psi_in_completo = kron(psi_in_completo, kron(psi_in_t, psi_in_s));
        case 2
            psi_in_completo = kron(psi_in_c, psi_in_t);
            psi_in_completo = kron(psi_in_completo, kron(psi_in_a, psi_in_s));
        case 3
            psi_in_completo = kron(psi_in_c, psi_in_t);
            psi_in_completo = kron(psi_in_completo, kron(psi_in_s, psi_in_a));
        otherwise
            disp('ERRORE:: indice del qudit da correggere sbagliato');
    end

%   Creiamo e vettoriziamo la matrice di densità iniziale del sistema
%   totale
    rho = vettorizza_matrice(psi_in_completo * psi_in_completo');
end