function [gl_completo, gl_ridotto] = genera_gate_logico(sim, ang)
%   Funzione che dati i parametri del gate logico da implementare ne
%   costruisce l'operatore relativo

    I = speye(sim.dim_a);
    gate_fisico = rotazione_piana(ang(1), ang(2));

%   Scriviamo il gate logico ridotto in base computazionale
    gl_ridotto  = sim.CB * kron(gate_fisico, speye(sim.dim_q/2)) * sim.CB';

%   Scriviamo il gate logico completo in base computazionale
    gl_completo  = kron(sim.CB,I) * kron(kron(gate_fisico, speye(sim.dim_q/2)),I) * kron(sim.CB',I);
end