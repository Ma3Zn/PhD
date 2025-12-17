function cw = pulisci_codewords(cw)
%   Funzione che date in input codewords ottenute numericamente forza la
%   disgiunzione dei supporti e rinormalizza i vettori

%   Recuperiamo le dimensioni del sistema
    dim = length(cw(:,1));

%   Recuperiamo gli indici degli elementi di modulo maggiore della
%   codewords
    [~,I] = sort(abs(cw(:,1)),'descend');

%   Settiamo a zero tutti gli elementi al di fuori del supporto delle
%   codewords
    cw(I(dim/2+1:end),1) = 0;
    
%   Normalizziamo la codewords
    cw(:,1) = cw(:,1) / norm(cw(:,1),2);

%   Ripetiamo per la seconda codewords. Recuperiamo gli indici degli
%   elementi di modulo maggiore della codewords
    [~,I] = sort(abs(cw(:,2)),'descend');

%   Settiamo a zero tutti gli elementi al di fuori del supporto delle
%   codewords
    cw(I(dim/2+1:end),2) = 0;
    
%   Normalizziamo la codewords
    cw(:,2) = cw(:,2) / norm(cw(:,2),2);
end