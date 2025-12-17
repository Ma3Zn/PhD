function [x_init] = genera_condizione_iniziale_MOD(dim)
%   Funzione che data la dimensione del sistema genera la condizione
%   iniziale per il CSP considerato

%   Condizione iniziale |0>, |1>
    x_0 = zeros(dim,1);
    x_1 = zeros(dim,1);

    x_0(1) = 1;
    x_1(2) = 1;

    x_init = [x_0;x_1];

% % % %   Condizione iniziale equal superosition
% % %     x = ones(dim,1)/sqrt(dim);
% % %     
% % %     x_init = [x;x];

% % % %   Condizione iniziale supporti divisi
% % %     x_0 = zeros(dim,1);
% % %     x_1 = zeros(dim,1);
% % % 
% % %     x_0(1:2:end) = 1;
% % %     x_1(2:2:end) = 1;
% % % 
% % %     x_0 = x_0 / sqrt(dim/2);
% % %     x_1 = x_1 / sqrt(dim/2);
% % % 
% % %     x_init = [x_0;x_1];

% % % %   Condizione iniziale supporti corretti
% % %     global ERR
% % % 
% % %     x_0 = zeros(dim,1);
% % %     x_1 = zeros(dim,1);
% % %     
% % %     x_0(ERR.supp0) = 1;
% % %     x_1(ERR.supp1 - dim) = 1;
% % % 
% % %     x_0 = x_0 / norm(x_0,2);
% % %     x_1 = x_1 / norm(x_1,2);
% % % 
% % %     x_init = [x_0;x_1];
end