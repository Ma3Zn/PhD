function err = calcola_errore_tomografia_logica_qudit(sim, rho)
%   Funzione che esegue una tomografia logica del quit
    
%   Ricriviamo la rho attuale in base logica
    rho = sim.CB' * full(rho) * sim.CB;

%   Misura lungo Z
    z_bar = trace(rho * sim.ZL);

%   Misura lungo X
    x_bar = trace(rho * sim.XL);

%   Misura lungo Y
    y_bar = trace(rho * sim.YL);

    rho_L = 1/2 * [1 + z_bar, x_bar - 1i*y_bar ; ...
                   x_bar + 1i*y_bar, 1 - z_bar];

%   Calcoliamo l'errore visto dalla misura logica del sistema
    err = 1 - sim.psi_fin_L' * rho_L * sim.psi_fin_L;
end