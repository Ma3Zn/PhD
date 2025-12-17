function [Hs, Hz, Hxy] = costruisci_hamiltoniane(N, J, periodic)
%   Funzione che genera i termini di interazione heisenberg tra due spin
%   consecutivi e l'Hamiltoniana da cui preparare il ground state.
%
%   In particolare:
%       -   Hz{i}  = Sz_{i} * Sz_{i+1}
%       -   Hxy{i} = Sx_{i} * Sx_{i+1} + Sy_{i}*Sy_{i+1}
%
%   Ho scelto questa strategia per avere già tutto impostato per utilizzare
%   dei J differenti tra i vari spin dell'anello.

%   Dimensione del sistema
    dim = 2^N;

%   Allochiamo memoria per Hs, Hz ed Hxy
    Hs  = sparse(dim, dim);
    Hz  = sparse(dim, dim);
    Hxy = sparse(dim, dim);

%   Generiamo gli operatori di spin
    Sx = sparse(sop(1/2,'x'));
    Sy = sparse(sop(1/2,'y'));
    Sz = sparse(sop(1/2,'z'));

%   Iteriamo sulla catena di spin e costruiamo le nostre Hamiltoniane
    for it = 1:N-1

%       Calcoliamo lo spazio a sinistra del primo spin
        dim_sx_i = 2^(it - 1);
%       Calcoliamo lo spazio a destra del primo spin
        dim_dx_i = 2^(N - it);

%       Calcoliamo lo spazio a sinistra del secondo spin
        dim_sx_ii = 2^(it);
%       Calcoliamo lo spazio a destra del primo spin
        dim_dx_ii = 2^(N - it - 1);

%       Creiamo le opportune matrici identità
        I_sx_i  = speye(dim_sx_i);
        I_dx_i  = speye(dim_dx_i);
        I_sx_ii = speye(dim_sx_ii);
        I_dx_ii = speye(dim_dx_ii);

%       Espandiamo gli operatori di spin alla base prodotto
        Sx_i  = kron(I_sx_i , kron(Sx, I_dx_i ));
        Sx_ii = kron(I_sx_ii, kron(Sx, I_dx_ii));

        Sy_i  = kron(I_sx_i , kron(Sy, I_dx_i ));
        Sy_ii = kron(I_sx_ii, kron(Sy, I_dx_ii));

        Sz_i  = kron(I_sx_i , kron(Sz, I_dx_i ));
        Sz_ii = kron(I_sx_ii, kron(Sz, I_dx_ii));

        hz = 0.4;
        hz = 0;

%       Inizializziamo Hs{it}
        Hs = Hs + hz * Sz_i;

%       Inizializziamo Hz{it}
        Hz = Hz + Sz_i * Sz_ii;

%       Inizializziamo Hxy{it}
        Hxy = Hxy + Sx_i * Sx_ii + Sy_i * Sy_ii;
    end
    
    if (periodic == true)
%       Calcoliamo lo spazio a sinistra del primo spin
        dim_sx_i = 2^(N-1);
%       Calcoliamo lo spazio a destra del primo spin
        dim_dx_i = 1;

%       Calcoliamo lo spazio a sinistra del secondo spin
        dim_sx_ii = 1;
%       Calcoliamo lo spazio a destra del primo spin
        dim_dx_ii = 2^(N - 1);
    
%       Creiamo le opportune matrici identità
        I_sx_i  = speye(dim_sx_i);
        I_dx_i  = speye(dim_dx_i);
        I_sx_ii = speye(dim_sx_ii);
        I_dx_ii = speye(dim_dx_ii);
    
%       Espandiamo gli operatori di spin alla base prodotto
        Sx_i  = kron(I_sx_i , kron(Sx, I_dx_i ));
        Sx_ii = kron(I_sx_ii, kron(Sx, I_dx_ii));
    
        Sy_i  = kron(I_sx_i , kron(Sy, I_dx_i ));
        Sy_ii = kron(I_sx_ii, kron(Sy, I_dx_ii));
    
        Sz_i  = kron(I_sx_i , kron(Sz, I_dx_i ));
        Sz_ii = kron(I_sx_ii, kron(Sz, I_dx_ii));
    
        hz = 0.2;
%         hz = 0;
    
%       Inizializziamo Hs{it}
        Hs = Hs + hz * Sz_i;

%       Inizializziamo Hz{it}
        Hz = Hz + Sz_i * Sz_ii;
    
%       Inizializziamo Hxy{it}
        Hxy = Hxy + Sx_i * Sx_ii + Sy_i * Sy_ii;
    end

%   Scaliamo i termini di heisenberg per la costante di accoppiamento J
    Hxy = J * Hxy;
    Hz  = J * Hz;
end