function p = genera_distribuzione_errore_misura(p_err_mis, dim, k)
%   Funzione che dato in input l'errore di misura, la dimensione del
%   sistema e l'indice del proiettore generalizzato da considerare crea la
%   distribuzione delle probabilità d'errore della misura, i.e., storce il
%   proiettore

%%  Distribuzione linearmente decrescente

% % % %   Numero proiettori sx
% % %     n_sx = k - 1;
% % % 
% % % %   Numero proiettori dx
% % %     n_dx = dim - k;
% % % 
% % % %   Generiamo un vettore contenente i pesi di ogni proiettore
% % %     pesi = zeros(1,dim);
% % % 
% % %     if (n_dx > n_sx)
% % % 
% % %         it = 1;
% % % 
% % %         for j = dim:-1:k+1
% % %             pesi(j) = it;
% % %             it = it + 1;
% % %         end
% % % 
% % %         pesi(n_sx:-1:1) = pesi(k+1:k+n_sx);
% % %     else
% % % 
% % %         for j = 1:n_sx
% % %             pesi(j) = j;
% % %         end
% % % 
% % %         pesi(k+1:end) = pesi(n_sx:-1:n_sx-n_dx+1);
% % %     end
% % % 
% % % %   Calcoliamo in quante parti dividere l'errore di misura
% % %     s = sum(pesi);
% % % 
% % % %   Calcoliamo l'errore unitario di misura
% % %     err = p_err_mis / s;
% % % 
% % % %   Generiamo la distribuzione di errori di misura opportuna partendo dal
% % % %   vettore dei pesi
% % %     p = pesi * err;
% % % 
% % % %   Aggiungiamo alla distribuzione la probabilità di eseguire una misura
% % % %   corretta
% % %     p(k) = 1 - p_err_mis;

%%  Distribuzione costante

%   inizializziamo il vettore dei pesi
    pesi = ones(1,dim);
    pesi(k) = 0;

%   Calcoliamo la somma dei pesi
    s = dim - 1;

%   Calcoliamo l'errore oportuno
    err = p_err_mis / s;

%   Assegnamo le opportune probabilità
    p = pesi * err;

%   Aggiungiamo alla distribuzione la probabilità di eseguire una misura
%   corretta
    p(k) = 1 - p_err_mis;

%% Sanity Check finale

%   Verifichiamo che la distribuzione generata sia una distribuzione di
%   probabilità
    if (abs(sum(p) - 1) > 1e2*eps)
        fprintf(stdout, "Distribuzione pesi proiettori generalizzati NON valida\n\n");
    end
end