function x_init = genera_condizione_iniziale_STD(dim)
%   Funzione che data in input la dimensione del sistema calcola il punto
%   iniziale dell'algoritmo. Per valutare in prima battuta le performance
%   dell'algoritmo di ottimizzazione useremo i supporti già individuati in
%   lavori precedenti con tecniche differenti. 
    
    s = dim/2;
    x_init = [];

    switch(s)
        case 2
            supp0 = [1 4];
            supp1 = [2 3];

        case 3
            supp0 = [1 3 6];
            supp1 = [2 4 5];

        case 4
            supp0 = [1 2 3 8];
            supp1 = [4 5 6 7];

        case 5
            supp0 = [2 3 5 7 10];
            supp1 = [1 4 6 8 9 ];

        case 6
            supp0 = [2 4 5 7  9 12];
            supp1 = [1 3 6 8 10 11];

        otherwise
            fprintf('\n\n\nbllllllllllllllllllllllllllllllll\n\n\n\n');
            return;
    end

    x_0 = zeros(dim,1);
    x_1 = zeros(dim,1);

    x_0(supp0) = 1 / sqrt(s);
    x_1(supp1) = 1 / sqrt(s);

    x_init = [x_0;x_1];
end