function sim = genera_op_QEC(dim, idx_q)
%   Funzione che genera tutti gli operatori indipendenti dei parametri
%   delle simulazioni. L'indice [idx_q] fornito in input indica a quale dei
%   tre qubit logici applicare un'effettiva procedura di EC e a quale di
%   questi applicare una stabilizzaione istantanea. Accorgimento necessario
%   questo, in quanto non ho le risorse al momento per eseguire una
%   simulazione con EC incoerente su tutti e tre gli oggetti
%   contemporaneamente. [Possibile alternativa implemntare una EC
%   incoerente consecutiva sui tre oggetti]

%   Ricaviamo le dimensioni degli oggetti componenti il sistema
    dim_q = dim;
    dim_a = dim / 2;

%   Settiamo precisione e tolleranza per il calcolo dell'esponenziale
    prec = 1e-8;
    tol  = 1e-14;

%   Recuperiamo la matrice Gamma
    gamma = leggi_gamma(dim_q);

%   Creiamo le codewords del codice (che saranno anche la nostra matrice di
%   cambio base [logica -> computazionale])
    CB = sparse(genera_codewords(dim_q));

%   Creiamo la matrice del cambiamento di base [logica -> computazionale]
%   per il sistema completo
    CB_completa = crea_base_completa(CB, idx_q, dim_a);

%   Creiamo la matrice del cambiamento di base [logica -> computazionale]
%   per il sistema ridotto [QQS]
    CB_ridotta = kron(kron(CB, CB), CB);

%   Calcoliamo la matrice di riordinamento del qudit
    P = calcola_matrice_riordinamento(dim_q, dim_a);

%   Operatore di evoluzione incoerente completo NON scalato di T/T2
    sop_inc_completo = costruisci_superoperatori_incoerenti_completi(gamma, dim_q, dim_a, idx_q);

%   Operatore di evoluzione incoerente ridotto NON scalato di T/T2
    sop_inc_ridotto = costruisci_superoperatori_incoerenti_ridotti(gamma, dim_q);

%   Generiamo, in base computazionale, il superoperatore del commutatore
%   dell'operazione di ED da implementare [CU x CU x CU]
    sop_comm_CU = genera_superoperatore_comm_CU(dim_q, CB_completa, idx_q);

%   Generiamo, in base computazionale, il superoperatore del commutatore
%   delle operazioni di recovery (sul sistema ridotto)
    [R, sop_comm_R] = genera_superoperatore_comm_R(dim_q, CB_ridotta);

%   Generiamo i proiettori per il sistema totale
    prj_tot = genera_proiettori_totali(dim_q, CB_completa, idx_q);

%   Generiamo i proiettori per il sisteam ridotto
    prj_rid = genera_proiettori_ridotti(dim_q, CB_ridotta);

%   Riempiamo sim con i parametri necessari
    sim.prec                = prec;             % precisione per esponenziale
    sim.tol                 = tol;              % tolleranza per esponenziale
    sim.dim_q               = dim_q;            % dimenqione qudit
    sim.dim_a               = dim_a;            % dimensione ancilla
    sim.gamma               = gamma;            % matrice gamma
    sim.CB                  = CB;               % matrice cambio base qudit
    sim.CB_ridotta          = CB_ridotta;       % matrice cambio base sistema ridotto
    sim.CB_completa         = CB_completa;      % matrice cambio base sistema completo
    sim.P                   = P;                % matrice per lo scambio di qudit e ancilla
    sim.sop_inc_completo    = sop_inc_completo; % parte incoerente superoperatore evoluzione sistema completo
    sim.sop_inc_ridotto     = sop_inc_ridotto;  % parte incoerente superoperatore evoluzione sistema ridotto
    sim.sop_comm_CU         = sop_comm_CU;      % parte coerente superoperatore evoluzione CU [sistema completo]
    sim.R                   = R;                % operatori di recovery [sistema ridotto]
    sim.sop_comm_R          = sop_comm_R;       % parte coerente superoperatore recovery [sistema ridotto]
    sim.prj_tot             = prj_tot;          % proiettori ancille sistema completo
    sim.prj_rid             = prj_rid;          % proiettori ideali sistema ridotto
end