function [H0, Hz, H_heis] = costruisci_hamiltoniane(spin_chain, periodic)
%   Funzione che genera le due Hamiltoniane estreme della nostra evoluzione
%   di partenza. In particolare H0 è l'Hamiltoniana di partenza con
%   interazione Ising mentre H1 è l'hamiltoninaa finale con interazione di
%   Heisenberg

%   Dimensione del sistema
    dim = 1;
    dim_vals = zeros(length(spin_chain), 1);

%   Numero di spin
    n = length(spin_chain);

%   Calcoliamo la dimensione del sistema
    for it = 1:n
        dim_vals(it) = 2*spin_chain(it) + 1;
        dim = dim * dim_vals(it);
    end

%   Allochiamo memoria per H0 ed H_z
    H0 = zeros(dim, dim);
    Hz = zeros(dim, dim);

%   Allochiamo lo spazio per tutti i termini che ci serviranno per il
%   calcolo dell'Hamiltoniana tempo dipendente
    H_heis.Sx_odd = zeros(dim, dim);
    H_heis.Sy_odd = zeros(dim, dim);
    H_heis.Sz_odd = zeros(dim, dim);

    H_heis.Sx_even = zeros(dim, dim);
    H_heis.Sy_even = zeros(dim, dim);
    H_heis.Sz_even = zeros(dim, dim);

% % % %   DBG
% % %     H_heis.i  = zeros(dim, dim);

%   Iteriamo sulla catena di spin e costruiamo le nostre Hamiltoniane
    for it = 1:n-1

%       Costruiamo gli operatori di spin nella base profotto
        Sx_i  = crea_operatore_prodotto(spin_chain, dim_vals, it,   'x');
        Sx_ii = crea_operatore_prodotto(spin_chain, dim_vals, it+1, 'x');

        Sy_i  = crea_operatore_prodotto(spin_chain, dim_vals, it,   'y');
        Sy_ii = crea_operatore_prodotto(spin_chain, dim_vals, it+1, 'y');

        Sz_i  = crea_operatore_prodotto(spin_chain, dim_vals, it,   'z');
        Sz_ii = crea_operatore_prodotto(spin_chain, dim_vals, it+1, 'z');

        hz = 0.04 * 0;

%       Aggiungiamo il termine corrente ad H0
        H0 = H0 + Sz_i * Sz_ii + hz * Sz_i;

%       Aggiungiamo il termine corrente ad H1
        Hz = Hz + hz * Sz_i;

%       Sommiamo i vari operatori per poi poter calcolare H_heis
        if (mod(it,2))
            H_heis.Sx_odd = H_heis.Sx_odd + Sx_i * Sx_ii;
            H_heis.Sy_odd = H_heis.Sy_odd + Sy_i * Sy_ii;
            H_heis.Sz_odd = H_heis.Sz_odd + Sz_i * Sz_ii;
        else
            H_heis.Sx_even = H_heis.Sx_even + Sx_i * Sx_ii;
            H_heis.Sy_even = H_heis.Sy_even + Sy_i * Sy_ii;
            H_heis.Sz_even = H_heis.Sz_even + Sz_i * Sz_ii;
        end
    end
    
%   Controlliamo se vogliamo una catena aperta o chiusa
    if (periodic && length(spin_chain) > 2)

%       Costruiamo gli operatori di spin nella base prodotto
        Sx_i  = crea_operatore_prodotto(spin_chain, dim_vals, n, 'x');
        Sx_ii = crea_operatore_prodotto(spin_chain, dim_vals, 1, 'x');

        Sy_i  = crea_operatore_prodotto(spin_chain, dim_vals, n, 'y');
        Sy_ii = crea_operatore_prodotto(spin_chain, dim_vals, 1, 'y');

        Sz_i  = crea_operatore_prodotto(spin_chain, dim_vals, n, 'z');
        Sz_ii = crea_operatore_prodotto(spin_chain, dim_vals, 1, 'z');

        hz = 0.6;

%       Aggiungiamo il termine corrente ad H0
        H0 = H0 + Sz_i * Sz_ii + hz * Sz_i;

%       Aggiungiamo il termine corrente ad Hz
        Hz = Hz + hz * Sz_i;

%       Sommiamo i vari operatori per poi poter calcolare H_heis
        if (mod(n,2))
            H_heis.Sx_odd = H_heis.Sx_odd + Sx_i * Sx_ii;
            H_heis.Sy_odd = H_heis.Sy_odd + Sy_i * Sy_ii;
            H_heis.Sz_odd = H_heis.Sz_odd + Sz_i * Sz_ii;
        else
            H_heis.Sx_even = H_heis.Sx_even + Sx_i * Sx_ii;
            H_heis.Sy_even = H_heis.Sy_even + Sy_i * Sy_ii;
            H_heis.Sz_even = H_heis.Sz_even + Sz_i * Sz_ii;
        end 
    end
end