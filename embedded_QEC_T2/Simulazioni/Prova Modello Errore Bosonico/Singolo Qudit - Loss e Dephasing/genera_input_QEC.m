function [sim] = genera_input_QEC(dim)
%   Funzione che genera tutti i dati di input necessari alla simulazione e
%   li salva all'interno di una opportuna struttura dati sim

%   Creiamo le codewords del codice (che saranno anche la nostra matrice di
%   cambio base [logica -> computazionale])
    CB = sparse(genera_codewords());

%   Creaiamo i proiettori del codice
    prj = genera_proiettori(dim);

%   Creiamo i proiettori sul solo spazio del qudit
    prj_q = genera_proiettori_qudit(dim);

%   Creiamo i proiettori logici per la misura logica del qubit
%   virtualizzato
    prj_l = genera_proiettori_logici_qudit(dim^2/2, dim);

%   Generiamo, in base computazionale, il superoperatore del commutatore
%   dell'operazione di ED da implementare [CU x CU x CU]
    [CU, sop_comm_CU] = genera_superoperatore_comm_CU(dim, kron(CB, speye(dim/2)));

%   Generiamo, in base computazionale, il superoperatore del commutatore
%   delle operazioni di recovery (sul sistema ridotto)
    [R, sop_comm_R] = genera_superoperatore_comm_R(dim, CB);

    for k = 1 : dim/2
        prj{k}   = kron(CB, eye(dim/2)) * prj{k} * kron(CB', eye(dim/2));
        prj_q{k} = CB * prj_q{k} * CB';

%       Rimuoviamo la sporcizia numerica
        prj{k}   = pulisci_matrice(  prj{k}, 1e3 * eps);
        prj_q{k} = pulisci_matrice(prj_q{k}, 1e3 * eps);
    end

    prj_l{1} = CB * prj_l{1} * CB';
    prj_l{2} = CB * prj_l{2} * CB';

%   Operatori necessari per tomografia logica del qudit  
    Z  = [1 0; 0 -1];
    ZL = kron(Z,eye(dim/2));

    X  = [0 1; 1 0];
    XL = kron(X,eye(dim/2));

    Y  = [0 -1i; 1i 0];
    YL = kron(Y,eye(dim/2));
    
    H = Hadamard(1/2);
    HL = kron(H, eye(dim/2));

    Rx = [1 -1i; -1i 1] / sqrt(2);
    RL = kron(Rx, eye(dim/2));

%   Riempiamo sim con i parametri necessari
    sim.dim_q   = dim;      % dimensione del qudit
    sim.dim_a   = dim/2;    % dimensione dell'ancilla
    sim.dim     = dim^2/2;  % dimensione del sistema
    sim.CB      = CB;       % matrice del cambio di base [log -> comp]
    sim.CU      = CU;       % stabilizzatore
    sim.prj     = prj;      % proiettori "k"
    sim.prj_q   = prj_q;    % proiettori ideali sul solo spazio del qudit
    sim.prj_l   = prj_l;    % proiettori logici del qudit
    sim.ZL      = ZL;       % Z logico
    sim.XL      = XL;       % X logico
    sim.YL      = YL;       % Y logico
    sim.HL      = HL;       % Gate logico per misure logiche in base X
    sim.RL      = RL;       % Gate logico per misure logiche in base Y
    sim.T2      = 0;        % T2
    sim.T_tot   = 0;        % tempo simulazione

%   Costruiamo gli operatori d'errore del sistema [in base computaizonale]
    E = cell(1,2);

    tmp = diag(sqrt(1:dim-1),1);
    E{1} = kron(tmp, eye(dim/2));

    tmp = tmp' * tmp;
    E{2} = kron(tmp/pi, eye(dim/2));

    tmp = diag(sqrt(1:dim/2-1),1);
    E{3} = kron(eye(dim), tmp);

    tmp = tmp' * tmp;
    E{4} = kron(eye(dim), tmp/pi);
    
    sim.E = E;

%   Operatore di evoluzione incoerente completo NON scalato di T/T2
    sop_inc_completo = costruisci_superoperatori_incoerenti(sim);

%   Operatore di evoluzione incoerente ridotto NON scalato di T/T2
    sop_inc_ridotto = costruisci_superoperatori_incoerenti_ridotti(sim);

    sim.sop_inc_completo    = sop_inc_completo; % parte incoerente superoperatore evoluzione sistema completo
    sim.sop_inc_ridotto     = sop_inc_ridotto;  % parte incoerente superoperatore evoluzione sistema ridotto

    sim.sop_comm_CU         = sop_comm_CU;      % parte coerente superoperatore evoluzione CU [sistema completo]
    sim.sop_comm_R          = sop_comm_R;       % parte coerente superoperatore evoluzione CU [sistema ridotto]
    sim.R                   = R;                % operatori di recovery [sistema ridotto]

% %
% % %   Generiamo le condizioni iniziali
% %     global alpha
% %     global beta
% % 
% %     psi_in_0              = zeros(dim, 1);
% %     psi_in_0(1)           = alpha;
% %     psi_in_0(dim/2 + 1)   = beta;
% % 
% %     psi_in_1 = zeros(dim/2,1);
% %     psi_in_1(1) = 1;
% % 
% %     psi_in = kron(psi_in_0, psi_in_1);
% %     psi_in = kron(CB, eye(dim/2)) * psi_in;
% % 
% % %   Condizioni iniziali del sistema
% %     rho = vettorizza_matrice(psi_in * psi_in');
% % 
% % %   Scriviamo il gate logico in base computazionale
% %     global gate_fisico
% %     gate_logico = kron(gate_fisico, eye(dim/2));
% % 
% % %   Generiamo le condizioni finali
% %     psi_fin = gate_logico * psi_in_0;
% %
% % %   Eseguaimo opportune trasformazioni per fare la misura in una base
% % %   logica arbitraria
% %     psi_fin = HL * psi_fin;
% %     psi_fin = RL * psi_fin;
% %
% % %   Scriviamo il gate logico e le condizioni finali in base computazionale
% %     psi_fin     = CB * psi_fin;
% % 
% %     sim.rho     = rho;      % condizioni iniziali [vettorizzata]
% %     sim.psi_fin = psi_fin;  % stato finale ESATTO
% %
% %     sim.gate_logico_ridotto = gate_logico;           % gate logico del sistema ridotto
% %     sim.gate_logico = kron(gate_logico, eye(dim/2)); % gate logico del sistema completo
% %
end