function [gl_c, gl_r] = genera_gate_logico_errato(sim, phi, q)
%   Funzione che dati in input i parametri della simulazione e l'angolo
%   del C-phi logico da implementare genera in base computazionale la forma
%   dell'operatore opportuno

    dim_c = sim.dim_a * (sim.dim_q/2)^3 * 2;
    dim_r = (sim.dim_q/2)^3 * 2;

%   Costruiamo la rappresentazione ridotta del nostro operatore, un CC-phi
%   (con target lo switch)
%     Rz  = rotazione_Z(phi);
    Rz = rotazione_piana(pi/2,phi);
    
    CCP = eye(4);

%   Generiamo il gate logico per il sistema completo scritto in base logica
    I_q = speye(sim.dim_q/2);
    I_q_a = speye(sim.dim_q/2 * sim.dim_a);
    switch q
        case 1
            gl_c = kron(I_q_a, kron(I_q, kron(CCP, I_q)));
            gl_c(end+1:end+dim_c, end+1:end+dim_c) = kron(I_q_a, kron(I_q, kron(Rz, I_q)));
        case 2
            gl_c = kron(I_q, kron(I_q_a, kron(CCP, I_q)));
            gl_c(end+1:end+dim_c, end+1:end+dim_c) = kron(I_q, kron(I_q_a, kron(Rz, I_q)));
        case 3
            gl_c = kron(I_q, kron(I_q, kron(CCP, I_q_a)));
            gl_c(end+1:end+dim_c, end+1:end+dim_c) = kron(I_q, kron(I_q, kron(Rz, I_q_a)));
        otherwise
            disp('ERRORE:: indice del qudit da correggere sbagliato');
    end

%   Rimuoviamo eventuale sporcizia numerica
    gl_c  = pulisci_matrice(gl_c, 1e3 * eps);

%   Generiamo il gate logico per il sistema ridotto scritto in base logica
    I_q = speye(sim.dim_q/2);
%     gl_r = kron(I_q, kron(I_q, kron(CCP, I_q)));
    gl_r = kron(I_q, kron(I_q, kron(CCP, I_q)));
    gl_r(end+1:end+dim_r, end+1:end+dim_r) = kron(I_q, kron(I_q, kron(Rz, I_q)));

%   Riscriviamo il gate logico in base computazionale
    gl_r = sim.CB_ridotta * gl_r * sim.CB_ridotta';

%   Rimuoviamo eventuale sporcizia numerica
    gl_r = pulisci_matrice(gl_r, 1e3 * eps);
end