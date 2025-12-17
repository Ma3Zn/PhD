function sop_err = costruisci_superoperatori_incoerenti_ridotti(sim)
%   Funzione che genera i superopertori incoerenti del sistema dati in
%   input gli opportuni parametri

%   Variabile d'appoggio
    z = sparse(sim.dim_q, sim.dim_q);

%   Superoperatore errore
    sop_err = sparse(sim.dim_q^2, sim.dim_q^2);

    for mu = 1:sim.dim_q
        for nu = 1:sim.dim_q

            A = z;
            B = z;

            A(mu,mu) = 1;
            B(nu,nu) = 1;

%           Calcoliamo il termine sempre presente
            tmp = 2 * genera_superoperatore_sx(A) ...
                    * genera_superoperatore_dx(B);

            if (mu == nu)
                tmp = tmp - genera_superoperatore_sx(A) ...
                          - genera_superoperatore_dx(B);
            end

%           Scaliamo il termine attuale e sommiamolo al superoperatore
%           totale
            sop_err = sop_err + sim.gamma(mu, nu) * tmp;
        end
    end
end