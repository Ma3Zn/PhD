function sop_err = costruisci_superoperatori_incoerenti_ridotti(gamma, dim_q)
%   Funzione che genera i superopertori incoerenti del sistema dati in
%   input gli opportuni parametri

    dim = dim_q;

%   Superoperatore errore
    sop_err = sparse(dim^6, dim^6);

%   Recuperiamo le dimensioni del sistema
    dims = [dim_q, dim_q, dim_q];

%   Cicliamo sui vari oggetti costituenti il sistema fisico
    for i = 1:3
        %   Variabile d'appoggio
        z = sparse(dims(i), dims(i));
%       Costruiamo il superoperatore d'errore
        for mu = 1:dims(i)
            for nu = 1:dims(i)
    
                A = z;
                B = z;
    
                A(mu,mu) = 1;
                B(nu,nu) = 1;

                I_sx = speye(prod(dims(1:i-1)));
                I_dx = speye(prod(dims(i+1:end)));
                
                A = kron(kron(I_sx, A), I_dx);
                B = kron(kron(I_sx, B), I_dx);
    
%               Calcoliamo il termine sempre presente
                tmp = 2 * genera_superoperatore_sx(A) ...
                        * genera_superoperatore_dx(B);
    
                if (mu == nu)
                    tmp = tmp - genera_superoperatore_sx(A) ...
                              - genera_superoperatore_dx(B);
                end
    
%               Scaliamo il termine attuale e sommiamolo al superoperatore
%               totale
                sop_err = sop_err + gamma(mu, nu) * tmp;
            end
        end
    end
end