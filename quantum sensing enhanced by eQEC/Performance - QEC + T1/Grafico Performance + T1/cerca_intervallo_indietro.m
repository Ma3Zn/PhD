function [Bx_a, Bx_b] = cerca_intervallo_indietro(S, Bx, T2, T_rot, lv_pop, T_max)
%   Funzione per la ricerca di un intervallo di valori di Bx in cui
%   eseguire una ricerca dicotomica. Il valore passato di Bx non è
%   misurabile

    mis = false;

%   Settiamo il valore iniziale delle variabili per la ricerca dicotomica
    Bx_b = 2*Bx;
    Bx_a = Bx;

%   Aggiorniamo il T_max
    T_max = T_max / 2;

%   Iniziamo la ricerca dicotomica
    while(~mis)
%       Verifichiamo la detectabilità del campo attuale
        [mis, T_mis] = verifica_detectabilita_campo(S, Bx_b, T2, T_rot, lv_pop, T_max);

%       Controlliamo se il protocollo è riuscito o meno a trovare
%       l'intersezione
        if (~mis)
%           Proviamo con un Bx più grande
            Bx_a = Bx_b;
            Bx_b = 2 * Bx_b;

%           Aggiorniamo il T_max stimato per il nuovo valore di Bx
            T_max = T_max / 2;
        end
    end
end