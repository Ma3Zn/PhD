function A = genera_operatori_base(W, dim)
%   Funzione test per la scrittura degli operatori base della tomografia
%   per un qudit di spin (dim-1)/2

%   Recuperiamo un valore di spin fittizio
    spin = (dim - 1)/2;

%   Allochiamo lo spazio per i nostri operatori
    A = cell(1, dim^2);

%   Teniamo traccia del numero di operatori creato
    idx = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Base -sperimentale-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %   Operatore da cui far partire la costruzione della base
%     E = zeros(dim,dim);
% 
% %   Riempiamo la diafonale inferiore di E
%     for i = 1:dim-1
%         E = E + diag(ones(dim - i,1), i);
%     end
% 
% %   Normalizziamo E
%     E = E / sqrt(trace(E'*E));
% 
% %   Generiamo la base ON
%     A = genera_base_sperimentale(E, dim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Stevens
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %   Operatori di Steevens
%     for k = 0:(dim - 1)
%         for q = -k:k
%             op_sv = stev(spin,k,q);
% 
%             A{idx} = op_sv / norm(full(op_sv), 'fro');
%             idx = idx + 1;
%         end
%     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   |n><m|
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for n = 1:dim
        for m = 1:dim
            A{idx} = zeros(dim,dim);
            A{idx}(n,m) = 1;
    
            idx = idx + 1;
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   |w><y|
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %   Usiamo come base degli operatori, gli operatori banali che si ottengono
% %   come ket-bra dagli autovettori di W
% 
% %   ATTENZIONE:: NON è una base ortonormale --> GS??
% 
% %   Diagonalizziamo W
%     [V, ~] = eig(W);
% 
% %   Costruiamo gli operatori
%     for n = 1:dim
%         for m = 1:dim
%             A{idx} = V(:,n)' * V(:,m);
%             A{idx} = A{idx} / norm(A{idx}, 'fro');
%             
%             idx = idx + 1;
%         end
%     end
end