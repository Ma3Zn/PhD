function sop_inc = calcola_sop_inc_T1_ridotto(sim)
%   Funzione per la generazione dei superoperatori della parte incoerente
%   della Lindblad master equation

%   ATTENZIONE:: |a> ed |b> in questa scrittura sono gli autostati
%   dell'Hamiltoniana. quindi questo lindbladiano è uguale sia in
%   interaction che in schrodinger picture

%   Recuperiamo la dimensione del qudit
    dim = sim.dim_q;

%   Variabile d'appoggio
    z = sparse(dim, dim);

%   Superoperatore errore
    sop_inc = sparse(dim^2, dim^2);

    for a = 1:dim
        for b = 1:dim

            AA = z;
            AB = z;
            BA = z;

            AA(a,a) = 1;
            AB(a,b) = 1;
            BA(b,a) = 1;
            
            tmp = genera_superoperatore_sx(BA) * genera_superoperatore_dx(AB) ...
                -(genera_superoperatore_sx(AA) + genera_superoperatore_dx(AA))/2;

%           Scaliamo il termine attuale e sommiamolo al superoperatore
%           totale
            sop_inc = sop_inc + sim.W(a,b) * tmp;
        end
    end
end