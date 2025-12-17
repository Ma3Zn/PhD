function [M, a] = genera_codewords()
%   Funzione che genera le codewords per la protezione di un photon loss
%   sulle note di quanto proposto da Girvin ma modificate per avere una
%   dimensione appropriata

    a = diag(sqrt(1:7),1);
    n = a'*a;

%%  Secondo me sono sbagliate!! Verifica e correggi

    M(:,1) = [1 0 0 0 sqrt(3) 0 0 0].'/2;
    M(:,2) = [0 0 0 1 0 0 0 0].';
    M(:,3) = [sqrt(3) 0 0 0 -1 0 0 0].'/2;
    M(:,3) = (sqrt(3)*M(:,1) - M(:,3))/2;
    M(:,3) = gs(M(:,1:2), M(:,3));

    M(:,4) = [0 0 sqrt(3) 0  0 0 1 0].'/2;
    M(:,5) = a * M(:,4);
    M(:,5) = M(:,5) / norm(M(:,5));
    M(:,6) = [0 0 1 0 0 0 -sqrt(3) 0].'/2;
    M(:,6) = (sqrt(3)*M(:,4) - M(:,6))/2;
    M(:,6) = gs(M(:,4:5), M(:,6));

    M(:,7) = gs(M,[1 1 1 1 1 0 0 0].');
    M(:,8) = [0 0 0 0 0 0 0 1].';

    MM = M;

    M(:,4) = MM(:,7);
    M(:,5) = MM(:,4);
    M(:,6) = MM(:,5);
    M(:,7) = MM(:,6);

%%  METTILE TUTTE:: TBD
    fprintf('Condizioni di Knill-Laflamme:\n');
    M(:,1)'*a'*a*M(:,1) - M(:,5)'*a'*a*M(:,5)
    M(:,1)'*a*M(:,1) - M(:,5)'*a*M(:,5)
    M(:,1)'*a'*M(:,1) - M(:,5)'*a'*M(:,5)

    M(:,1)'*a'*a*M(:,5)
    M(:,1)'*a'*M(:,5)
    M(:,1)'*a*M(:,5)
end