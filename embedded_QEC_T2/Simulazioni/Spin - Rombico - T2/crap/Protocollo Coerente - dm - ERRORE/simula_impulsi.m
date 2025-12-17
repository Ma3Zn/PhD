function [sim, freq] = simula_impulsi(sim)
%   funzione che approssimia numericamente l'evoluzione del sistema durante
%   gli impulsi oscillanti lungo Z e calcola la frequenza di oscillazione
%   logica

%   Rendiamo disponibile alla funzione per il calcolo del forzante tutti i
%   dati della simulazione da eseguire
    global param;
    param = sim;

%   Dimensione del sistema
    dim = sim.sys.dim;

%   Mesh temporale
    t_fin   = 4e6; %ns
    Ntime   = 5e4;
    time    = linspace(0, t_fin, Ntime);
    
%   Salviamo la mesh temporale utilizzata
    sim.time = time;

%   Allochiamo lo spazio per l'evoluzione della matrice di densità
    sim.evol   = zeros(dim, dim, Ntime+1);
    sim.evol_L = sim.evol;

%   Scriviamo lo stato iniziale del sistema in base computazionale
    sim.evol_L(:,:,1) = sim.sys.psi_in * sim.sys.psi_in';

    psi_in = sim.sys.CB * sim.sys.psi_in;
    sim.evol(:,:,1)   = psi_in * psi_in';

% % % %   DBG: Scriviamo lo stato iniziale in base computazionale
% % %     sim.evol(:,:,1)   = sim.sys.psi_in * sim.sys.psi_in';
% % %     sim.evol_L(:,:,1) = sim.sys.CB' * sim.evol(:,:,1) * sim.sys.CB;

    %   Matrice di densità iniziale del sistema
    rho0 = sim.evol(:,:,1);

%   Reshape della condizine iniziale
    y0 = reshape(rho0, dim^2, 1);

%   Evoluzione del sistema
    tic
    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
    [t,y] = ode45(@forzante_lindblad,time,y0,options);
    toc

%   Reshape delle soluzioni approssimate
    for i = 1:Ntime
%       Costruiamo l'operatore di cambio di base tra interaction picture e
%       rotating frame
%         U = genera_operatore_rf(sim.imp, t(i));

        sim.evol(:,:,i+1)   = reshape(y(i,:).',dim,dim);
%         sim.evol(:,:,i+1)   = U * reshape(y(i,:).',dim,dim) * U';
        sim.evol_L(:,:,i+1) = sim.sys.CB' * sim.evol(:,:,i) * sim.sys.CB;
    end

%   Conduziamo un analisi di fourier per ottenere la frequenza di
%   oscillazione logica
    
    % TBD
    freq = 0;
end