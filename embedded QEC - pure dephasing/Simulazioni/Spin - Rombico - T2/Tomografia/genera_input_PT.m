function input = genera_input_PT(dim)
%   Funzione che genera gli input necessari alla process tomography del CU
%   nell'algoritmo di QEC.

%   Leggiamo i parametri per il T1
    W = leggi_W(dim);

%   Leggiamo i parematri per il T2
    g = leggi_g(dim);

%   Durata della fase di idle
    t_U = 5e1;

%   T1 del sistema
    T1 = 5e6;

%   T2 del sistema
    T2 = 1e5;

%   Generiamo gli operatori su cui fare la tomografia
    A = genera_operatori_base(W, dim);

%   Numero di oggetti che compongono il sistema
    n_corpi = 1;

%   Gate di cui fare la tomografia
    U = speye(dim);

%   Generiamo il superoperatore incoerente per il T1
    sop_inc_T1 = costruisci_sop_inc_T1(W, dim);

%   Generiamo il superoperatore incoerente per il T2
    sop_inc_T2 = costruisci_sop_inc_T2(g, dim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Parametri richiesti nella variabile di output di questa funzione per
%   l'esecuzione della tomografia della procedura proposta
    input.n_corpi = n_corpi;        % Numero di oggetti che compongono il sistema
    input.dim = dim;                % Dimensione del sistema
    input.t   = t_U;                % Tempo di esecuzione del gate soggetto a tomografia
    input.T1  = T1;                 % Tempi di decadimento degli oggetti del sistema
    input.T2  = T2;                 % Tempi di coerenza degli oggetti del sistema
    input.U   = U;                  % Gate da sottoporre a tomografia
    input.A   = A;                  % Operatori su cui eseguire la tomografia
    input.sop_inc_T1 = sop_inc_T1;  % Superoperatore incoerente del T1
    input.sop_inc_T2 = sop_inc_T2;  % Superoperatore incoerente del T2
end