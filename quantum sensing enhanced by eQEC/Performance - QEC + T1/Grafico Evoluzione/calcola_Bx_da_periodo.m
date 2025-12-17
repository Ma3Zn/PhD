function Bx_mis = calcola_Bx_da_periodo(sim)
%   Dato il semi-periodo di rabi recuperiamo il valore del campo trasverso
%   nell'ipotesi in cui la QEC abbia ridotto l'errore a sufficienza
%   affinche sin e cos si incrocino ancora nel semiperiodo

%   T = 2*K(2,1) / (Bz(2,1) * sim.Bx * (sim.g * sim.mu_b)^2 / (sim.h_bar^2))
%   K(i,j) = sim.w(i,j) * gen_U(i,j) / comm_zx(i,j)

%   Recuperiamo gli indici delle transizioni da eseguire
    r = find(sim.CB(:,              1))';
    c = find(sim.CB(:,sim.dim_q/2 + 1))';

%   Calcoliamo Bx dal semiperiodo per ogni transizione eseguita
    n = 1;
    for i = r
        for j = c
            Bx(n) = 2*sim.K(i,j) / ( sim.T * sim.Bz(i,j) * (sim.g * sim.mu_b)^2 / (sim.h_bar^2) );
            n = n + 1;
        end
    end

%   Mediamo le stime di Bx
    Bx_mis = somma_Kahan(Bx) / (n-1);
end