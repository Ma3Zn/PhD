function [H, S] = costruisci_matrici_SbE_ZNE(exp_val)
% Funzoine che passati i valori di aspettazione di N potenze
% dell'Hamiltoniana ne costruisce le relative matrici di una Krylov-like
% SbE

%   Recuperiamo il numero di potenze dell'Hamiltoniana da usare
    N = length(exp_val);

%   Calcoliamo la dimensione della SbE [il +1 è perchè usiamo anche I] 
    dim = (N - 1)/2 + 1;

%   Allochiamo le matrici della SbE
    H = zeros(dim);
    S = zeros(dim);

%   Forziamo a mano una soluzione che sappiamo funzoinare per {I,H}
    H(1,1) = exp_val(1);
    H(1,2) = exp_val(2);
    H(2,1) = exp_val(2);
    H(2,2) = exp_val(3);

    S(1,1) = 1;
    S(1,2) = exp_val(1);
    S(2,1) = exp_val(1);
    S(2,2) = exp_val(2);

%   DA SISTEMARE
% % % %   Costruiamo le matrici con gli opportuni valori di aspettazione
% % %     for r = 1:dim
% % %         for c = 1:dim
% % %             idx = (r - 1) * dim + c;
% % % 
% % %             H(r,c) = exp_val(idx);
% % % 
% % %             if (r==1 && c==1)
% % %                 S(r,c) = 1;
% % %             else
% % %                 S(r,c) = exp_val(idx-1);
% % %             end
% % %         end
% % %     end
end