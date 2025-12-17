function rho_fin = esegui_QEC(sim)
%   Funzione che applica un ciclo di QEC al sistema

    rho_fin = zeros(sim.dim_q^2,1);

%   Estendiamo lo stato del sistema per includere anche l'ancilla
    rho = vettorizza_matrice(kron(devettorizza_matrice(sim.rho), sim.rho_a));

%   Applichiamo il CU al sistema
    rho = sim.sop_CU * rho;

% % % %   Calcoliamo le probabilità di selezionare un k
% % %     for j = 1:sim.dim_q/2
% % %         sim.p(j) = trace( sim.prj{j} * rho );
% % %     end

%   Eseguiamo una fase di IDLE per simulare il processo di misura
    rho = devettorizza_matrice(sim.sop_idle * rho);

%   Portiamoci dietro una miscela delle rho
    for k = 1:sim.dim_a
%       Proiettiamo lo stato del sistema
        rho_k = sim.prj{k} * rho * sim.prj{k}';
        
%       Tracciamo fuori l'ancilla
        rho_k = vettorizza_matrice(partial_trace_1(rho_k, sim.dim_q, sim.dim_a));

%       Eseguiamo lo step di recovery
        rho_fin = rho_fin + sim.sop_R{k} * rho_k;
    end
end