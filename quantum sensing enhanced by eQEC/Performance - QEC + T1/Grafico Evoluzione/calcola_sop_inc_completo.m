function sop_err = calcola_sop_inc_completo(sim)
%   Funzione che genera i superopertori incoerenti del sistema dati in
%   input gli opportuni parametri

%   Superoperatore errore
    sop_err = sparse(sim.dim^2, sim.dim^2);

%   Recuperiamo le dimensioni del sistema
    dims = [sim.dim_q, sim.dim_a];

%   Cicliamo sui vari oggetti costituenti il sistema fisico
    for i = 1:2
        %   Variabile d'appoggio
        z = sparse(dims(i),dims(i));
%       Costruiamo il superoperatore d'errore
        for mu = 1:dims(i)
            for nu = 1:dims(i)
    
                A = z;
                B = z;
    
                A(mu,mu) = 1;
                B(nu,nu) = 1;

                if i == 1
                    A = kron(A, speye(sim.dim_a));
                    B = kron(B, speye(sim.dim_a));
                else
                    A = kron(speye(sim.dim_q), A);
                    B = kron(speye(sim.dim_q), B);
                end
    
%               Calcoliamo il termine sempre presente
                tmp = 2 * genera_superoperatore_sx(A) ...
                        * genera_superoperatore_dx(B);
    
                if (mu == nu)
                    tmp = tmp - genera_superoperatore_sx(A) ...
                              - genera_superoperatore_dx(B);
                end
    
%               Scaliamo il termine attuale e sommiamolo al superoperatore
%               totale
                sop_err = sop_err + sim.gamma(mu, nu) * tmp;
            end
        end
    end
end