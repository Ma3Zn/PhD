function [sim, sop] = calcola_sop_comm_rot(sim)
%   Funzione che genera il superoperatore del commutatore del generatore
%   della rotazione logica
%
%   Questo generatore è quello dato (in ns) nel sistema di riferimento
%   rotante opportuno mandando impulsi polarizzati circolarmente

%   Allocazione variabili
    comm_gen = zeros(sim.dim_q, sim.dim_q);
    Bz = zeros(sim.dim_q, sim.dim_q);
    K  = zeros(sim.dim_q, sim.dim_q);

%   Recuperiamo gli indici delle transizioni da eseguire
    r = find(sim.CB(:,              1))';
    c = find(sim.CB(:,sim.dim_q/2 + 1))';

%   Calcoliamo il logaritmo dell'unitaria da implementare
    gen_U = pulisci_matrice(1i*logm(sim.GL), 1e3 * eps);
    gen_U = sim.CB * gen_U * sim.CB';

%   Calcoliamo il commutatore tra Mx ed Mz diagonale
    comm_zx = diag(diag(sim.Mz))*sim.Mx - sim.Mx*diag(diag(sim.Mz));

%   Calcoliamo i vari K
    for i = r
        for j = c
            K(i,j) = sim.w(i,j) * gen_U(i,j) / comm_zx(i,j);
        end
    end

%   Salviamo K all'interno di sim
    sim.K = K;

%   Recuperiamo il K più grande --> transizione con elemento di matrice più
%   piccolo
    [M, R] = max(K);
    [~, C] = max(M);

    idx_K_max = [R(C), C];

%   Calcoliamo i vari valori di Bz (in questo modo avremo un drive
%   corretto (contemporaneo) dell'operazione logica da eseguire)
    for i = r
        for j = c
            if (i ~= idx_K_max(1) || j ~= idx_K_max(2))
                Bz(i,j) = K(i,j) / K(idx_K_max(1),idx_K_max(2)) * sim.Bz_max;
            else
                Bz(i,j) = sim.Bz_max;
            end
        end
    end

%   Simmetrizziamo la matrice dei Bz ( che non ha alcun elemento diagonale)
    Bz = Bz + Bz';

%   Salviamo Bz all'interno di sim
    sim.Bz = Bz;

%   Calcoliamo gli elementi di matrice del nostro generatore
    for i = 1:sim.dim_q
        for j = 1:sim.dim_q
            if i ~= j
                comm_gen(i,j) = Bz(i,j) * sim.Bx ...
                              * ((sim.g * sim.mu_b)^2 / (sim.h_bar^2 * 2*sim.w(i,j))) ...
                              * comm_zx(i,j);
            end
        end
    end

%   Scaliamo la durata temporale consistentemente con il numero di
%   operatori di cui dover fare la trotterizzazione -- 
%   2S+1/2  <=> 1/(sim.dim_q/2) <=> 1/sim.dim_a
    comm_gen = comm_gen / sim.dim_a;

%   T = 2 * K(2,1) / (Bz(2,1) * sim.Bx * (sim.g * sim.mu_b)^2 / (sim.h_bar^2))
%   T*comm_gen = gen_U
    sim.T_est = sim.dim_a * 2 * K(2,1) / (Bz(2,1) * sim.Bx * (sim.g * sim.mu_b)^2 / (sim.h_bar^2));

%   Creiamo il superoperatore commutatore
    sop = genera_superoperatore_sx(comm_gen) ...
        - genera_superoperatore_dx(comm_gen) ;

%   Moltiplichiamo opportunamente per l'unità immaginaria
    sop = -1i * sop;
end