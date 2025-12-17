function sop_inc = costruisci_sop_inc_T2(g, dim)
%   Funzione per la generazione dei superoperatori della parte incoerente
%   della Lindblad master equation

%   ATTENZIONE:: |a> e |b> in questa scrittura sono gli autostati
%   dell'Hamiltoniana. quindi questo lindbladiano è uguale sia in
%   interaction che in schrodinger picture

%   Variabile d'appoggio
    z = sparse(dim, dim);

%   Superoperatore errore
    sop_inc = sparse(dim^2, dim^2);

    for a = 1:dim
        for b = 1:dim

            AA = z;
            BB = z;

            AA(a,a) = 1;
            BB(b,b) = 1;

            tmp = genera_superoperatore_sx(AA) * genera_superoperatore_dx(BB);

            sop_inc = sop_inc + g(a, b) * tmp;
        end
    end

    sop_inc = sparse(sop_inc);
end