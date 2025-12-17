function [Bx, sens, T_max] = cerca_in_avanti(S, Bx, sens, step, T2, T_rot, lv_pop, T_max)
%   Funzione per la ricerca di un intervallo di valori di Bx in cui
%   eseguire una ricerca dicotomica. Il valore passato di Bx è misurabile
    
%   Variabili d'appoggio
    Bx_old = 0;
    sens_old = inf;
    T_max_old = inf;

%   Contatore
    iter = 1;

%   Iniziamo la ricerca dicotomica
    while(sens < sens_old)

        Bx_old_old = Bx_old;
        sens_old_old = sens_old;
        T_max_old_old = T_max_old;
        
        Bx_old = Bx;
        sens_old = sens;
        T_max_old = T_max;

%       Settiamo il nuovo campo da testare
        Bx = (1 - step) * Bx_old;

%       Aggiorniamo il T_max
        T_max = (1 + step) * T_max_old;

%       Verifichiamo la detectabilità del campo attuale
        [~, T_mis, sens] = verifica_detectabilita_campo(S, Bx, T2, T_rot, lv_pop, T_max);

%       Aggiorniamo il contatore
        iter = iter + 1;
    end

%   Abbiamo superato il minimo ritorniamo il valore opportuno dell-ultimo
%   campo testato prima del minimo
    if (iter > 2)
        Bx = Bx_old_old;
        sens = sens_old_old;
        T_max = T_max_old_old;
    else
        Bx = Bx_old;
        sens = sens_old;
        T_max = T_max_old;
    end
end