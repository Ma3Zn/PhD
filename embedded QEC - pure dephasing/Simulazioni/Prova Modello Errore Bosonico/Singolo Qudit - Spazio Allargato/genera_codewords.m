function [M, a] = genera_codewords()
%   Funzione che genera le codewords per la protezione di un photon loss
%   sulle note di quanto proposto da Girvin ma modificate per avere una
%   dimensione appropriata

    M(:,1) = [0 0 1 0 0 0 0 0].';
    M(:,2) = [0 1 0 0 0 0 0 0].';
    M(:,3) = [1 0 0 0 1 0 0 0].'/sqrt(2);
    M(:,4) = [0 0 0 1 0 0 0 0].';
    M(:,5) = gs(M,[1 1 1 1 0 0 0 0].');
    M(:,6) = [0 0 0 0 0 1 0 0].';
    M(:,7) = [0 0 0 0 0 0 1 0].';
    M(:,8) = [0 0 0 0 0 0 0 1].';

%   Riordiniamo opportunamente i vettori della base
    MM = M;
    M(:,3) = MM(:,5);
    M(:,4) = MM(:,7);
    M(:,5) = MM(:,3);
    M(:,6) = MM(:,4);
    M(:,7) = MM(:,6);


    a = diag(sqrt(1:7),1);

    fprintf('Condizioni di Knill-Laflamme:\n');
    M(:,1)'*a'*a*M(:,1) - M(:,5)'*a'*a*M(:,5)
    M(:,1)'*a*M(:,1) - M(:,5)'*a*M(:,5)
    M(:,1)'*a'*M(:,1) - M(:,5)'*a'*M(:,5)

    M(:,1)'*a'*a*M(:,5)
    M(:,1)'*a'*M(:,5)
    M(:,1)'*a*M(:,5)
end