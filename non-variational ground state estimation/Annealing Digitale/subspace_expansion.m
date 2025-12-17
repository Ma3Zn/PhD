function gs = subspace_expansion(sim, psi, n_potenze)
%       Funzione che dati i parametri della simulazione, gli stati iniziali
%       a cui applicare le potenze delle Hamiltoiane e l'esponente massi a
%       cui elevare tali operatori genera il sottospazio di krylov
%       relativo, risolve il problema agli autovettori generalizzato come
%       in "Quantum Power Method by a Superposition of Time-Evolved States"
%       e ritorna l'approssimazione del ground state calcolata come
%       sovrapposizione degli stati base del sottospazio

%   Allochiamo la memoria per il ground state approssimato
    gs = zeros(length(psi(:,1)),1);

%   Costruiamo l'Hamiltoniana di Heisenberg del sistema target
    J = sim.J;
    H = J.z*sim.H_heis.Sz_odd  + J.x*sim.H_heis.Sx_odd  + J.y*sim.H_heis.Sy_odd ...
      + J.z*sim.H_heis.Sz_even + J.x*sim.H_heis.Sx_even + J.y*sim.H_heis.Sy_even;

%   Calcoliamo la dimensione effettiva dello spazio di Krylov
    dim = n_potenze + 1;

%   Allochiamo lo spazio per i vettori dello spazio di Krylov
    K_vector = zeros(length(psi(:,1)), dim);

%   Allochiamo lo spazio opportuno per le matrici del problema
%   generalizzato agli autovalori all'interno dello spazio di Krilov
    H_krylov = zeros(dim);
    S_krylov = zeros(dim);

%   Generiamo i vettori dello spazio di Krylov
    for j = 0:n_potenze
%       Calcoliamo la potenza dell'Hamiltoniana da applicare
        Hj = H^j;

%       Generiamo l'opportuno stato di Krylov
        K_vector(:,j+1) = Hj * psi;
    end

%   Costruiamo le matrici H ed S del prblema agli autovettori
%   generalizzato attuale
    for k = 1:dim
        for p = 1:dim
            H_krylov(k,p) = K_vector(:,k)' * H * K_vector(:,p);
            S_krylov(k,p) = K_vector(:,k)' * K_vector(:,p);
        end
    end

%   Risolviamo il problema agli autovalori generalizzato Hv = eSv con
%   le routine di MATLAB, dettagli su come risolvere tale problema
%   possono essere trovati in un qualsiasi testo di analisi numerica o
%   anche nell'articolo "Quantum Power Method by a Superposition of
%   Time-Evolved States"
    [V, D] = eig(real(H_krylov),real(S_krylov));

%   Recuperiamo l'autovettore relativo al minimo autovalore del
%   problema generalizzato
    [~, idx] = min(diag(D));
    v_min = V(:,idx);

%   Ricosruiamo l'approssimazione del ground state 
    for j = 1:dim
        gs = gs + v_min(j) * K_vector(:,j);
    end

%   Normalizziamo lo stato finale
    gs = gs / norm(gs);
end