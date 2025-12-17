function op = genera_O21(spin)
%   Funzione che genera loperatore SxSz + SzSx per un Giant spin con
%   interazione rombica

    h_bar = 6.582119569000000e-04;
    mu_b  = 0.057883826350000e+00;
    g = 2;

%%  Operatori di spin

    dim = 2*spin + 1;
    
    Sx = sop(spin, 'x');
    Sy = sop(spin, 'y');
    Sz = sop(spin, 'z');

%%  Parametri del sistema

    d = -0.1;
    e_max = d / 3;
    e = -0.03 /100;
    B = 0.35; % T
    a = 0;
    
    Ham = d * Sz * Sz + e * (Sx * Sx - Sy * Sy) + g * mu_b * B * (Sz * cos(a) + Sx * sin(a));
    Ham = Ham / (spin * (2*spin - 1));
    [V, H0] = eig(Ham);

%%  Scrittura di GAMMA

    Sxr = diag(V' * Sx * V);
    Szr = diag(V' * Sz * V);
    
    G = zeros(dim, dim);
    
    for n = 1:dim
        for m = 1:dim
            G(m,n) = Sxr(m)*Sxr(n) + Sxr(m)*Szr(n) + Szr(m)*Sxr(n) + Szr(m)*Szr(n);
        end
    end
    
    Sxr = V' * Sx * V;
    Syr = V' * Sy * V;
    Szr = V' * Sz * V;

%% Scrittura di W

    op = Sxr*Szr + Szr*Sxr;
    
    op = triu(abs(op).^2);
    op = op / norm(op, 'fro');

end