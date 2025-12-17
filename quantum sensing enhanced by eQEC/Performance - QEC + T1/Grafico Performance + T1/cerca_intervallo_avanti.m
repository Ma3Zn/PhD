function [Bx_a, Bx_b, T_max] = cerca_intervallo_avanti(S, Bx, T2, T_rot, lv_pop, T_max)
%   Funzione per la ricerca di un intervallo di valori di Bx in cui
%   eseguire una ricerca dicotomica. Il valore passato di Bx è misurabile

    mis = true;
    
%   Settiamo il valore iniziale delle variabili per la ricerca dicotomica
    Bx_b = Bx;
    Bx_a = Bx/2;

%   Aggiorniamo il T_max
    T_max = 2 * T_max;

%   Iniziamo la ricerca dicotomica
    while(mis)
%       Verifichiamo la detectabilità del campo attuale
        [mis, T_mis] = verifica_detectabilita_campo(S, Bx_a, T2, T_rot, lv_pop, T_max);

%       Controlliamo se il protocollo è riuscito o meno a trovare
%       l'intersezione
        if (mis)
%           Proviamo con un Bx più piccolo
            Bx_b = Bx_a;
            Bx_a = Bx_a / 2;

%           Aggiorniamo il T_max stimato per il nuovo valore di Bx
            T_max = 2 * T_mis;
        end
    end
end