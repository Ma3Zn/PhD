function [v] = random_superposition(n)
%   Funzione he genera una sovrapposizione arbitraria di stati    

%   Distribuzione randomica delle fasi
    angle = randn(n,1);

%   Distribuzione randomica delle popolazioni
    pop = abs(randn(n,1));

%   Normalizziamo le popolazioni
    pop = pop / norm(pop,'fro');
    
    v = pop .* exp(1i * angle);
end
