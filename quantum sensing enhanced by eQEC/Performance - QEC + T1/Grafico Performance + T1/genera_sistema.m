function sys = genera_sistema(S, Bx, T1, T2)
%   Funzione che genera una struttura dati contenente tutti i paramentri
%   del sistema corrente

%%   Definizione costanti
    h_bar = 6.582119569000000e-04;
    mu_b  = 0.057883826350000e+00;
    g = 2;

%%  Operatori di spin

    spin = S;
    dim = 2*spin + 1;
    
    Sx = sop(spin, 'x');
    Sy = sop(spin, 'y');
    Sz = sop(spin, 'z');

%%  Parametri hamiltoniana

    d = -0.1;
    e_max = d / 3;
    e = -0.03;
    B = 0.35; % T
    a = 0;
    
    Ham = d * Sz * Sz + e * (Sx * Sx - Sy * Sy) + g * mu_b * B * (Sz * cos(a) + Sx * sin(a));

    Ham = Ham / (S * (2*S - 1)); % Scaling opportuno di D ed E con S

    [V, H0] = eig(Ham);

    H0 = H0(end:-1:1, end:-1:1);
    V  =  V(   :    , end:-1:1);

%   Creiamo la matrice delle frequenze risonanti
    w = zeros(dim,dim);
    for i = 1:dim
        for j = 1:dim
            w(i,j) = (H0(i,i) - H0(j,j))/h_bar;
        end
    end

%%  Parametri campo oscillante
    Bz_max = 0.01; %T

%%  Momenti magnetici

    Mxr = (V' * Sx * V);
    Myr = (V' * Sy * V);
    Mzr = (V' * Sz * V);

%%  Path ai parametri corretti

    radice = strcat('./par/', num2str(dim), '/');

%%  Codewords
%   Apriamo gli stream su file

    nome_0 = strcat(radice, 'cw_0.bin');
    nome_1 = strcat(radice, 'cw_1.bin');

    fileID_0 = fopen(nome_0, 'r');
    fileID_1 = fopen(nome_1, 'r');

%   Leggiamo le due codifiche
    ew0 = fread(fileID_0, 'double');
    ew1 = fread(fileID_1, 'double');

%   Formattiamo correttamente le due codifiche
    M  = [reshape(ew0, dim, dim/2), reshape(ew1, dim, dim/2)];

%   Chiudiamo gli stream su file
    fclose(fileID_0);
    fclose(fileID_1);

%%  Recuperiamo gli indici delle popolazioni di |0L> ed |1L>
    
    supp0_v = calcola_indici_vettorizzati(find(M(:,         1))', dim);
    supp1_v = calcola_indici_vettorizzati(find(M(:, dim/2 + 1))', dim);

%%  Lettura Gamma per T2
    nome = strcat(radice, 'G.txt');

%   Leggiamo la matrice gamma
    G = load(nome);

%%  Lettura W per T1
    nome = strcat(radice, 'Fs.txt');

%   Leggiamo il fattore di scala
    Fs = load(nome);

    nome = strcat(radice, 'W.txt');

%   Leggiamo la matrice W
    W = load(nome);

%   Scaliamo opportunamente W
    W = Fs * W;

%%  Condizione iniziale qudit

%   inizializziamo il sistema in |0_L> [ base computazionale ]
    psi_in    = zeros(dim,1);
    psi_in(1) = 1;
    
    psi_in = M * psi_in;

    rho = vettorizza_matrice(psi_in * psi_in');

%%  Condizione iniziale ancilla
    psi_in_a    = zeros(dim/2,1);
    psi_in_a(1) = 1;

    rho_a = psi_in_a * psi_in_a';

%%  Generazione proiettori stabilizzazione
    prj = genera_proiettori(M, dim);

%%  Probabilità di stabilizzazione
    p = zeros(dim/2,1);

%%  Rotazione logica

%   Rappresentazione in base computazionale del gate logico da eseguire
    GL = kron(rotazione_piana(pi/2, pi/2), eye(dim/2));

%%  Tempi esecuzioni rotazioni logiche
    t_es(4) = 2.17;
    t_es(6) = 6.91;
    t_es(8) = 22.2;

%%  Creazione struttura dati

    sys.g       = g;       % costante
    sys.mu_b    = mu_b;    % costante
    sys.h_bar   = h_bar;   % costante
    sys.dim_q   = dim;     % dimensinoe del qudit
    sys.dim_a   = dim/2;   % dimensione dell'ancilla
    sys.dim     = dim^2/2; % dimensione totale del sistema [qudit+ancilla]
    sys.Mx      = Mxr;     % momento magnetico x
    sys.My      = Myr;     % momento magnetico y
    sys.Mz      = Mzr;     % momento magnetico z
    sys.Bz_max  = Bz_max;  % Massimo valore di Bz
    sys.Bx      = Bx;      % Valore del campo esterno da misurare
    sys.H0      = H0;      % Hamiltonaina del sistema
    sys.w       = w;       % matrice delle frequenze risonanti
    sys.CB      = M;       % matrice del cambio di base [comp <- log]
    sys.rho     = rho;     % stato iniziale del sistema [vettorizzato]
    sys.rho_old = rho;     % stato iniziale del sistema [vettorizzato]
    sys.rho_a   = rho_a;   % stato iniziale dell'ancilla
    sys.W       = W;       % matrice W per il T1
    sys.gamma   = G;       % matrice gamma per il T2
    sys.T1      = T1;      % T1 del sistema
    sys.T2      = T2;      % T2 del sistema
    sys.supp0_v = supp0_v; % indici vettorizzati popolazioni |0L>
    sys.supp1_v = supp1_v; % indici vettorizzati popolazioni |1L>
    sys.prj     = prj;     % proiettori di stabilizzazione
    sys.p       = p;       % probbailità di stabilizzazione
    sys.GL      = GL;      % rotazione logica in base computazionale
    sys.t_es    = t_es;    % tempo esecuzione operazioni logiche
    
end