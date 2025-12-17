function [flag, T_mis, sens] = verifica_detectabilita_campo(S, Bx, T2, T_rot, lv_pop, T_max)
%   Funzione che valuta se il campo passato è detectabile o meno

%   Valore di default del flag di riuscita o meno del protocollo
    flag = true;

%   Generiamo i parametri della simulazione
    sim = genera_parametri_simulazione(S, Bx, T2, T_rot);

%   Eseguiamo il protocollo di sensing
    [T_mis, T_es] = esegui_protocollo_sensing(sim, lv_pop, T_max);

%   Controlliamo se il protocollo è riuscito o meno a trovare
%   l'intersezione
    if (~(T_mis > 0))
        flag = false;
        sens = inf;
        return;
    end

%   Offset relativo per approssimazione derivata
    delta = 1e-10;

%   Se il campo è detectabile calcoliamo anche la sensibilità relativa
    Bx_sens = Bx + Bx*delta;

%   Generiamo i parametri della simulazione
    sim = genera_parametri_simulazione(S, Bx_sens, T2, T_rot);

%   Eseguiamo il protocollo di sensing
    T_tmp = esegui_protocollo_sensing(sim, lv_pop, T_max);

%   Approssimazione derivata
    dSdB = abs( (T_mis - T_tmp) / (Bx - Bx_sens) );
    derivata = 1/dSdB;

%   Calcolo approssimato della sensitività
    sens = derivata / sqrt(floor(1/(T_es*1e-9)));
end