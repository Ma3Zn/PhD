function CU = genera_CU(dim)
%   Funzione che fornito in input lo spin del qudit genera il gate CU
%   opportuno

%   Costruiamo CU in base logica
    CU = eye(dim/2);
    for k = 2:dim/2
        tmp = eye(dim/2);
        
        tmp(1,1) = 0;
        tmp(k,k) = 0;

        tmp(k,1) =  1;
        tmp(1,k) = -1;

        CU(end+1:end+dim/2,end+1:end+dim/2) = tmp;
    end

%   La rappresentazione di CU è indipendente da l
    CU(end+1:2*end,end+1:2*end) = CU;
end