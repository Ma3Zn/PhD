function Bx = ricerca_dicotomica_Bx_min(S, Bx_a, Bx_b, T2, T_rot, lv_pop, T_max, N_step)
%   Funzione che fornito in input un intervallo di valori in cui ricercare
%   il valor minimo del campo misurabile dal protocollo ricerca tale valore
%   tramite un algoritmo dicotomico e convergente

    mis = false;
    T_mis = T_max;

%   Stimiamo il campo minimo
    for i = 1:N_step

%       Calcoliamo il campo attuale da testare
        Bx = (Bx_a + Bx_b) / 2;

%       Aggiorniamo il valore di T_max in modo opportuno
        if (mis)
%           L'ultimo Bx era misurabile
%             T_max = 2*(1+log10(Bx/Bx_b)) * T_mis;
            T_max = T_mis * ( 1 / (Bx_b/Bx) );
        else
%           L'ultimo Bx NON era misurabile
%             T_max =   (1+log10(Bx_a/Bx)) * T_max;
            T_max = T_max * ( 1 / (Bx/Bx_a) );
        end

%       Verifichiamo la misurabilità del campo attuale
        [mis, T_mis] = verifica_detectabilita_campo(S, Bx, T2, T_rot, lv_pop, T_max);

%       Aggiorniamo gli estremi dell'intervallo
        if (mis)
%           Il campo attuale è misurabile
            Bx_b = Bx;
        else
%           Il campo attuale NON è misurabile
            Bx_a = Bx;
        end

    end
end