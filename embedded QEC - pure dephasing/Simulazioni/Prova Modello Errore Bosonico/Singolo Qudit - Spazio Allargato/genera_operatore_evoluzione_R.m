function evol_R = genera_operatore_evoluzione_R(sim, T, T2)
%   Funzione che genera il k-esimo superoperatore di evoluzione incoerente 
%   del sistema ridotto [Q Q S] che implementa l'opportuno step di recovery

%   Allochiamo lo spazio per evol_R
    evol_R = cell(1,sim.dim_a);

%   Scorriamo i vari k
    for k = 1:sim.dim_a

%       Aggiungiamo la parte incoerente
        comm_R = sim.sop_comm_R{k};
        comm_R = comm_R + sim.sop_inc_ridotto * (T / T2);

%       Esponenziamo l'operaore
        evol_R{k} = fastExpm(comm_R);
    end
end