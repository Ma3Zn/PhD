function [R] = genera_operatori_recovery(dim)
%   Funzione che dato uno spin ed il k in cui effetturare la recovery
%   genera gli opportuni operatori

    R = cell(1,dim/2);

    R{1} = eye(dim);

    for j = 2 : dim/2
        R{j} = eye(dim);

        R{j}(1,1) = 0;
        R{j}(j,j) = 0;

        R{j}(1,j) = 1;
        R{j}(j,1) = -1;

        R{j}(1+dim/2,1+dim/2) = 0;
        R{j}(j+dim/2,j+dim/2) = 0;

        R{j}(1+dim/2,j+dim/2) = 1;
        R{j}(j+dim/2,1+dim/2) = -1;
    end 
end