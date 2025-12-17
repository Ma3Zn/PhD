function dev = grafico_errore_medio(fig, dim, T2_vals, angle_vals, init_superposition)
%   Funzione che dati in input la dimensione del sistema, i valori di T2
%   per cui eseguire le simulazioni, i parametri delle rotazioni e gli
%   stati iniziali su cui mediare, calcola l'errore medio dovuto
%   all'esecuzione dei vari gate logici

%   Numero di T2 da simulare
    [~, n_T2] = size(T2_vals);

%   Numero di operazioni da testare
    [~, n_op] = size(angle_vals);

%   Numero di stati di input da provare
    [~, n_stati_iniziali] = size(init_superposition);

%   Errore medio della procedura
    errore_medio = zeros(3, n_T2);

%   Errore massimo della procedura
    errore_massimo = zeros(1, n_T2);

%   Errore minimo della procedura
    errore_minimo = ones(1, n_T2);

%   Cicliamo il qudit da correggere realmente    
    for q_corretto = 1:3
        
        fprintf('Posizione ancilla: %d\n', q_corretto);
        
        tic
        fprintf('Generazione operatori input\n');

    %   Generiamo l'input per la simulazione
        sim = genera_op_QEC(dim, q_corretto);
        toc

%       Indice di avanzamento (T2)
        idx_T2 = 1;

    %   Cicliamo sui valori di T2
        for T2 = T2_vals
            fprintf('\n1/T2: ');
            disp(1/T2);
            fprintf('\n');
            
    %       Inizializiamo l'indice di stato
            idx_op = 1;
    
    %       Inizializiamo l'errore di stato da mediare
            errore_op = zeros(1,n_op);
    
            tic
            fprintf('Generazione superoperatore evoluzione CU\n');

    %       Generiamo il superoperatore di evoluzione incoerente del CU
            comm_CU = sim.sop_comm_CU;
            comm_CU = comm_CU + sim.sop_inc_completo * (90/T2);
            sim.sop_evol_CU = fastExpm(comm_CU, sim.prec, sim.tol);
            
            toc

            tic
            fprintf('Generazione superoperatori evoluzione R\n');

    %       Generiamo il superoperatore di evoluzione incoerente delle recovery
            sim.sop_evol_R = genera_operatore_evoluzione_R(sim, 90, T2);
    
            toc

    %       Cicliamo sui parametri delle operazioni da testare
            for ang = angle_vals
                fprintf('\nphi: %f\n\n', ang);

    %           Scriviamo il gate logico in base computazionale
                [gl_completo, gl_ridotto]     = genera_gate_logico(sim, 0, q_corretto);
                
                tic
                fprintf('Generazione superoperatori evoluzione Gl\n');
    %           Calcoliamo il superoperatore sparso del gate logico (pulisci
    %           matrice torna già unamatrice sparsa)
                log_G  = sim.CB_completa * logm(full(gl_completo)) * sim.CB_completa';
                log_G  = pulisci_matrice(log_G, 1e3 * eps);
                comm_G = genera_superoperatore_sx(log_G) - genera_superoperatore_dx(log_G);
    
    %           Calcoliamo il superoperatore di evoluzioine
                comm_G = comm_G + sim.sop_inc_completo * (90/T2);
                sim.sop_evol_gl = fastExpm(comm_G, sim.prec, sim.tol);
    
                toc

    %           Inizializiamo l'indice di stato
                idx_stato = 1;
    
    %           Inizializiamo l'errore di stato da mediare
                errore_stato = zeros(1,n_stati_iniziali^2);
    
    %           Cicliamo sugli stati iniziali del control
                for amp_c = init_superposition
                    for amp_t = init_superposition
                        fprintf('\ncondizioni iniziali (C): ');
                        disp([amp_c(1) amp_c(2)]);
                        fprintf('condizioni iniziali (T)');
                        disp([amp_t(1), amp_t(2)]);
                        fprintf('\n');
    
    %                   Calcoliamo le condizioni iniziali del sistema
                        [sim.rho, psi_in] = crea_condizioni_iniziali(sim, amp_c, amp_t, q_corretto);
    
    %                   Calcoliamo le condizioni finali
                        sim.psi_fin = gl_ridotto * gl_ridotto * psi_in;

                        tic
                        fprintf('INIZIO procedura EC\n');
                        
    %                   Calcoliamo l'errore del circuito implementante l'algoritmo di QEC
                        [p_osservazione, errore] = calcola_errore(sim, q_corretto);
                        
                        fprintf('FINE procedura EC\n');
                        toc

    %                   Calcoliamo l'errore medio
                        errore_stato(idx_stato) = p_osservazione.' * errore;
    
    %                   Facciamo avanzare l'indice di stato
                        idx_stato = idx_stato + 1;
                    end
                end

    %           Recuperiamo il massimo ed il minimo valore per l'errore delle
    %           operazioni
                errore_massimo(idx_T2) = max(abs([errore_stato errore_massimo(idx_T2)]));
                errore_minimo(idx_T2)  = min(abs([errore_stato errore_minimo(idx_T2)]));
    
    %           Mediamo l'errore di stato
                errore_op(idx_op) = sum(errore_stato) / n_stati_iniziali^2;
    
    %           Facciamo avanzare l'indice di operazione
                idx_op = idx_op + 1;
            end
    
    %       Calcoliamo l'errore medio
            errore_medio(q_corretto, idx_T2) = sum(errore_op) / n_op;
    
    %       Facciamo avanzare l'indice di T2
            idx_T2 = idx_T2 + 1;
        end
    end

%   Mediamo gli errori medi riscontrati
    err_fin = (errore_medio(1,:) + errore_medio(2,:) + errore_medio(3,:)) / 3;

%   Facciamo il grafico della curva ottenuta
    fai_grafico(err_fin, T2_vals, 1e2 * eps);
    saveas(fig, './images/pseudo_thresholds.fig');

%   L'errore minimo ha poco senso
    dev = errore_massimo;
end
