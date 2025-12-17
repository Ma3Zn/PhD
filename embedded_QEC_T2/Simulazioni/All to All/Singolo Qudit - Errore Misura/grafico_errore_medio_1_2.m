function dev = grafico_errore_medio_1_2(T2_vals, angle_vals, init_superposition)
%   Funzione che dati in input i valori degli angoli delle rotazioni piane
%   da eseguire, i vari valori di T2 e gli stati iniziali su cui mediare la
%   fidelity fa il grafico della fidelity media con tanto di intorni
%   tubolari definiti dal valore massimo osservato

%   Numero di T2 da simulare
    [~, n_T2] = size(T2_vals);

%   Numero di operazioni da testare
    [~, n_op] = size(angle_vals);

%   Numero di stati di input da provare
    [~, n_stati_iniziali] = size(init_superposition);

%   Indice di avanzamento (T2)
    idx_T2 = 1;

%   Errore medio della procedura
    errore_medio = zeros(1, n_T2);

%   Errore massimo della procedura
    errore_massimo = zeros(1, n_T2);

%   Errore minimo della procedura
    errore_minimo = ones(1, n_T2);

%   Cicliamo sui valori di T2
    for T2 = T2_vals

%       Inizializiamo l'indice di stato
        idx_op = 1;

%       Inizializiamo l'errore di stato da mediare
        errore_op = zeros(1,n_op);

%       Cicliamo sui parametri delle operazioni da testare
        for R_par = angle_vals

%           Inizializiamo l'indice di stato
            idx_stato = 1;

%           Inizializiamo l'errore di stato da mediare
            errore_stato = zeros(1,n_stati_iniziali);

%           Cicliamo sugli stati iniziali del sistema
            for sup_in = init_superposition

%               Generiamo l'input per la simulazione (supponiamo la durata
%               del gate fisico di 10 ns)
                in = genera_input_1_2_med(T2, R_par, sup_in, 2);

%               Eseguiamo la simulazione
                rho_fin = devettorizza_matrice(in.sop_evol * in.rho);

%               Calcoliamo l'errore
                errore_stato(idx_stato) = 1 - in.psi_fin' * rho_fin * in.psi_fin;

%               Facciamo avanzare l'indice di stato
                idx_stato = idx_stato + 1;
            end
%           Recuperiamo il massimo ed il minimo valore per l'errore delle
%           operazioni
            errore_massimo(idx_T2) = max(abs([errore_stato errore_massimo(idx_T2)]));
            errore_minimo(idx_T2)  = min(abs([errore_stato errore_minimo(idx_T2)]));

%           Mediamo l'errore di stato
            errore_op(idx_op) = sum(errore_stato) / n_stati_iniziali;

%           Facciamo avanzare l'indice di operazione
            idx_op = idx_op + 1;
        end

%       Calcoliamo l'errore medio
        errore_medio(idx_T2) = sum(errore_op) / n_op;

%       Facciamo avanzare l'indice di T2
        idx_T2 = idx_T2 + 1;
    end

%   Facciamo il grafico della curva ottenuta
    fai_grafico(errore_medio.', T2_vals, 1e2 * eps);

% %   Ritorniamo gli intorni tubolari della curva ottenura
%     dev = [errore_massimo; errore_minimo];

    dev = errore_massimo;

end