function [flag, T_mis] = verifica_detectabilita_campo(S, Bx, T2, T_rot, lv_pop, T_max)
%   Funzione che valuta se il campo passato è detectabile o meno

%   Valore di default del flag di riuscita o meno del protocollo
    flag = true;

%   Generiamo i parametri della simulazione
    sim = genera_parametri_simulazione(S, Bx, T2, T_rot);

%   Eseguiamo il protocollo di sensing
    T_mis = esegui_protocollo_sensing(sim, lv_pop, T_max);

%   Controlliamo se il protocollo è riuscito o meno a trovare
%   l'intersezione
    if (~(T_mis > 0))
        flag = false;
    else
        T_mis = T_max;
    end
end