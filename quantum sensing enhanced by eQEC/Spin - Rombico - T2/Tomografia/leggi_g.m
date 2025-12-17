function g = leggi_g(dim)
%   Funzione per la lettura dei parametri per il T2

%   Leggiamo i ratei di T2
    G = load("par/G.txt");

%   Riduciamo le dimensioni di G
    G = G(1:dim, 1:dim);
    g = G;

%   Ricaviamo gamma da GAMMA
    for mu=1:dim
        for nu=1:dim
            g(mu,nu) = 2*G(mu,nu)-G(mu,mu)-G(nu,nu);
        end
    end
end