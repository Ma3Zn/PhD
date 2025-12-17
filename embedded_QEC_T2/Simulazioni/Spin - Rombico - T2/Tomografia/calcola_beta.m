function beta = calcola_beta(j, in)
%   Funzione che dati in inut la rho attuale e i parametri del sistema
%   calcola i rispettivi beta e li ritorna in una matrice di in.dim^2 righe
%   e in.dim^4 colonne.
%   j è l'indice di iterazione, necessario per la parallelizzazione del loop

%   Recuperiamo la dimensione del sistema
    dim = in.dim;

%   Recuperiamo gli indici m ed n di j
    m_j = rem(j, dim);

    if m_j == 0
        m_j = dim;
    end

    n_j = (j - m_j) / dim + 1;

%   Allochiamo beta
    beta = zeros(in.dim^2,in.dim^4);

%   Cicliamo sui vari k
    for n_k = 1:in.dim
        for m_k = 1:in.dim
%           Recuperiamo l'indice iterativo attuale
            k = (n_k-1) * in.dim + m_k;

%           Cicliamo sui vari operatori della base scelta
            for n = 1:in.dim^2

%               Costruiamo l'operatore relativo
                En = in.A{n};

                for m = 1:in.dim^2
%                   Costruiamo l'operatore relativo
                    Em = in.A{m}';

%                   Recuperiamo l'indice di riga di beta
                    riga = k;

%                   Recuperiamo l'indice di colonna di beta
                    colonna = (n-1) * in.dim^2 + m;

% %                   Calcoliamo beta_{(jk)(nm)} [ATTNEZIONE:: il ' su Em è
% %                   già stato messo quando inizializzato come in.A(:,:,m)']
%                     beta(riga,colonna) = trace(En*rho_j*Em*rho_k');

%                   Data la peculiare struttura di rho_k e rho_j, la
%                   traccia desiderata è ...
                    beta(riga, colonna) = En(n_k,n_j) * Em(m_j,m_k);
                end
            end
        end
    end

%   Ristrutturiamo la forma dell'outpt per eseguire la parallelizzazione del ciclo
    tmp = zeros(in.dim^4);
    tmp((j-1) * in.dim^2 + 1 : j * in.dim^2, :) = beta;
    beta = tmp;
end
