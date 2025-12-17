function sop_inc = calcola_sop_inc_T1_completo(sim)
%   Funzione per la generazione dei superoperatori della parte incoerente
%   della Lindblad master equation

%   ATTENZIONE:: |a> ed |b> in questa scrittura sono gli autostati
%   dell'Hamiltoniana. quindi questo lindbladiano è uguale sia in
%   interaction che in schrodinger picture

%   Superoperatore errore
    sop_inc = sparse(sim.dim^2, sim.dim^2);

%   Recuperiamo le dimensioni del sistema
    dims = [sim.dim_q, sim.dim_a];

%   Cicliamo sui vari oggetti costituenti il sistema fisico
    for i = 1:2

%       Variabile d'appoggio
        z = sparse(dims(i),dims(i));

%       Costruiamo il superoperatore d'errore
        for a = 1:dims(i)
            for b = 1:dims(i)

                AA = z;
                AB = z;
                BA = z;

                AA(a,a) = 1;
                AB(a,b) = 1;
                BA(b,a) = 1;

                if i == 1
                    AA = kron(AA, speye(sim.dim_a));
                    AB = kron(AB, speye(sim.dim_a));
                    BA = kron(BA, speye(sim.dim_a));
                else
                    AA = kron(speye(sim.dim_q), AA);
                    AB = kron(speye(sim.dim_q), AB);
                    BA = kron(speye(sim.dim_q), BA);
                end
            
                tmp = genera_superoperatore_sx(BA) * genera_superoperatore_dx(AB) ...
                    -(genera_superoperatore_sx(AA) + genera_superoperatore_dx(AA))/2;

%               Scaliamo il termine attuale e sommiamolo al superoperatore
%               totale
                sop_inc = sop_inc + sim.W(a,b) * tmp;
            end
        end
    end
end