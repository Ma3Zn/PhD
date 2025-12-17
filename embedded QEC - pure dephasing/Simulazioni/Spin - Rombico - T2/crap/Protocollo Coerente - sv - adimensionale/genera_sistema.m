function sys = genera_sistema()
%   Funzione che genera una struttura dati contenente tutti i paramentri
%   del sistema corrente

%%   Definizione costanti
    h_bar = 6.582119569e-04;
    mu_b  = 0.057883826350000e+00;

%%  Operatori di spin

    spin = 3/2;
    dim = 2*spin + 1;
    
    Sx = sop(spin, 'x');
    Sy = sop(spin, 'y');
    Sz = sop(spin, 'z');

%%  Parametri hamiltoniana

    d = 0.3;
    e_max = d / 3;
    e = e_max;
    g = 0.5; % g*mu_b * B0 => B0 ~= 5T;
    a = 0;
    
    Ham = d * Sz * Sz + e * (Sx * Sx - Sy * Sy) + g * (Sz * cos(a) + Sx * sin(a));
    [V, H0] = eig(Ham);

    H0 = H0(end:-1:1, end:-1:1);
    V  =  V(   :    , end:-1:1);

%   Scaliamo opportunamente gli autovalori di H0
    H0 = diag((diag(H0) - H0(1,1)));

%   Creiamo la matrice delle frequenze risonanti
    w = zeros(dim,dim);
    for i = 1:dim
        for j = 1:dim
            w(i,j) = H0(i,i) - H0(j,j);
        end
    end

%%  Parametri campo oscillante
    gz = 2;
    Bz_max = 0.02; %T

%%  Parametri campo statico incognito
    gx = 2;

%%  Momenti magnetici

    Mxr = (V' * Sx * V);
    Myr = (V' * Sy * V);
    Mzr = (V' * Sz * V);

%%  Codewords
%   Apriamo gli stream su file
    fileID_0 = fopen('./data/cw_4lv_0.bin', 'r');
    fileID_1 = fopen('./data/cw_4lv_1.bin', 'r');

%   Leggiamo le due codifiche
    ew0 = fread(fileID_0, 'double');
    ew1 = fread(fileID_1, 'double');

%   Formattiamo correttamente le due codifiche
    M  = [reshape(ew0, dim, dim/2), reshape(ew1, dim, dim/2)];

%   Chiudiamo gli stream su file
    fclose(fileID_0);
    fclose(fileID_1);

%%  Condizione iniziale sistema

%   inizializziamo il sistema in |0_L>
    psi_in = [1 0 0 0].';

%   DBG
%     psi_in = [4*sqrt(3) 2*sqrt(3)*1i -sqrt(3)*1i 1].'/8;
%     psi_in = [1 1 1 1].'/2;
%     psi_in = [0 0 1 0].';

%%  Rotazione logica

%   Rappresentazione in base logica del gate logico da eseguire
    GL = kron(rotazione_piana(pi/2, 0), eye(dim/2));

%%  Creazione struttura dati

    sys.mu_b   = mu_b;
    sys.h_bar  = h_bar;

    sys.dim     = dim;     % dimensine del sistema
    sys.Mx      = Mxr;     % momento magnetico x
    sys.My      = Myr;     % momento magnetico y
    sys.Mz      = Mzr;     % momento magnetico z
    sys.H0      = H0;      % Hamiltonaina del sistema
    sys.w       = w;       % matrice delle frequenze risonanti
    sys.CB      = M;       % matrice del cambio di base [comp <- log]
    sys.GL      = GL;      % gate logico da eseguire
    sys.psi_in  = psi_in;  % stato iniziale del sistema
    sys.gz      = gz;      % Accoppiamento lungo z
    sys.Bz_max  = Bz_max;  % Massimo valore di Bz
    sys.gx      = gx;      % Accoppiamento lungo x

end