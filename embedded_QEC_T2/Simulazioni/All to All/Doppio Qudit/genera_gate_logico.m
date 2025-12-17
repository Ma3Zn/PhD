function [gl_c, gl_r] = genera_gate_logico(sim, phi, q)
%   Funzione che dati in input i parametri della simulazione e l'angolo
%   del C-phi logico da implementare genera in base computazionale la forma
%   dell'operatore opportuno

    P = sim.P;

%%  Attenzione:: FUNZIONA SOLO PER dim=4
    
    R = rotazione_piana(pi/2, phi);
    
%     I0 = speye(eye(2));
%     I1 = speye(sim.dim_q/2);
% 
%     z1 = sparse(2*sim.dim_q/2,2*sim.dim_q/2);
%     z2 = kron(z1,z1);

%     gl_r_mod = kron([kron(I2, I2) kron(z2,z2); kron(z2,z2) [kron(I0, I1) z1; z1 kron(R, I1)]]);
    z = zeros(sim.dim_q);
    I = eye(sim.dim_q);
    R = kron(R, eye(sim.dim_q/2));
    
    if (sim.dim_q == 4)    
        gl_r = [I z z z z z z z z z z z z z z z;
                z I z z z z z z z z z z z z z z;
                z z I z z z z z z z z z z z z z;
                z z z I z z z z z z z z z z z z;
                z z z z I z z z z z z z z z z z;
                z z z z z I z z z z z z z z z z;
                z z z z z z I z z z z z z z z z;
                z z z z z z z I z z z z z z z z;
                z z z z z z z z I z z z z z z z;
                z z z z z z z z z I z z z z z z;
                z z z z z z z z z z R z z z z z;
                z z z z z z z z z z z R z z z z;
                z z z z z z z z z z z z I z z z;
                z z z z z z z z z z z z z I z z;
                z z z z z z z z z z z z z z R z;
                z z z z z z z z z z z z z z z R];

    elseif (sim.dim_q == 6)
        gl_r = [I z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z;
                z I z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z;
                z z I z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z;
                z z z I z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z;
                z z z z I z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z;
                z z z z z I z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z;
                z z z z z z I z z z z z z z z z z z z z z z z z z z z z z z z z z z z z;
                z z z z z z z I z z z z z z z z z z z z z z z z z z z z z z z z z z z z;
                z z z z z z z z I z z z z z z z z z z z z z z z z z z z z z z z z z z z;
                z z z z z z z z z I z z z z z z z z z z z z z z z z z z z z z z z z z z;
                z z z z z z z z z z I z z z z z z z z z z z z z z z z z z z z z z z z z;
                z z z z z z z z z z z I z z z z z z z z z z z z z z z z z z z z z z z z;
                z z z z z z z z z z z z I z z z z z z z z z z z z z z z z z z z z z z z;
                z z z z z z z z z z z z z I z z z z z z z z z z z z z z z z z z z z z z;
                z z z z z z z z z z z z z z I z z z z z z z z z z z z z z z z z z z z z;
                z z z z z z z z z z z z z z z I z z z z z z z z z z z z z z z z z z z z;
                z z z z z z z z z z z z z z z z I z z z z z z z z z z z z z z z z z z z;
                z z z z z z z z z z z z z z z z z I z z z z z z z z z z z z z z z z z z;
                z z z z z z z z z z z z z z z z z z I z z z z z z z z z z z z z z z z z;
                z z z z z z z z z z z z z z z z z z z I z z z z z z z z z z z z z z z z;
                z z z z z z z z z z z z z z z z z z z z I z z z z z z z z z z z z z z z;
                z z z z z z z z z z z z z z z z z z z z z R z z z z z z z z z z z z z z;
                z z z z z z z z z z z z z z z z z z z z z z R z z z z z z z z z z z z z;
                z z z z z z z z z z z z z z z z z z z z z z z R z z z z z z z z z z z z;
                z z z z z z z z z z z z z z z z z z z z z z z z I z z z z z z z z z z z;
                z z z z z z z z z z z z z z z z z z z z z z z z z I z z z z z z z z z z;
                z z z z z z z z z z z z z z z z z z z z z z z z z z I z z z z z z z z z;
                z z z z z z z z z z z z z z z z z z z z z z z z z z z R z z z z z z z z;
                z z z z z z z z z z z z z z z z z z z z z z z z z z z z R z z z z z z z;
                z z z z z z z z z z z z z z z z z z z z z z z z z z z z z R z z z z z z;
                z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z I z z z z z;
                z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z I z z z z;
                z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z I z z z;
                z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z R z z;
                z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z R z;
                z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z z R];
    end

%   Sparsifichiamo l'operatore del gate ottenuto
    gl_r = sparse(gl_r);

    I = speye(sim.dim_q^2);

%   Generiamo il gate logico completo (verificare correttezza)
    switch q
        case 1
            gl_c = kron(speye(sim.dim_a), gl_r);

%           Calcoliamo l'opportuno riordinamento della base
            riord = kron(P, I);
            
%           Riordiniamo la base scambiando ancilla e control [AQQS] --> [QAQS]
            gl_c = riord * gl_c * riord';
        case 2
            gl_c = kron(gl_r, speye(sim.dim_a));

%           Calcoliamo l'opportuno riordinamento della base
            riord = kron(I, P);

%           Riordiniamo la base scambiando ancilla e control [QQSA] --> [QQAS]
            gl_c = riord' * gl_c * riord;
        case 3
            gl_c = kron(gl_r, speye(sim.dim_a));
        otherwise
            disp('ERRORE:: indice del qudit da correggere sbagliato');
    end
    
%   Riscriviamo il gate logico in base computazionale
    gl_r = sim.CB_ridotta * gl_r * sim.CB_ridotta';

%   Rimuoviamo eventuale sporcizia numerica
    gl_r = pulisci_matrice(gl_r, 1e3 * eps);

end