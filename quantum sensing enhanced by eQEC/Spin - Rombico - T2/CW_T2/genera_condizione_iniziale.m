function x_init = genera_condizione_iniziale(dim)
%   Funzione che data in input la dimensione del sistema calcola il punto
%   iniziale dell'algoritmo. Per valutare in prima battuta le performance
%   dell'algoritmo di ottimizzazione useremo i supporti già individuati in
%   lavori precedenti con tecniche differenti. 

    x_0 = zeros(dim/2,1);
    x_1 = zeros(dim/2,1);

% % %     x_0(1) = 1;
% % %     x_1(1) = 1;

    x_0 = x_0 + 1 / sqrt(dim);
    x_1 = x_1 + 1 / sqrt(dim);

    x_init = [x_0; x_1];
end