function gamma = leggi_gamma(dim)
%   Funzione che data la dimensione dello spazio in input restituiscela
%   gamma opportuna

%   Leggiamo la matrice gamma
    nome = strcat('./par/', num2str(dim),'/G.txt');
    G = load(nome);

    gamma = G * 1e0;
end