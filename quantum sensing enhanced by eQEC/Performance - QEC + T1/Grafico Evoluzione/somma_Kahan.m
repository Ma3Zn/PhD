function somma = somma_Kahan(x)
%   Funzione per la somma cumulativa degli eleementi del vettore di input x
%   implementanta tramite l'algoritmo di Kaan-Babuska per la minimizzazione
%   dell'errore numerico

%   Variabile che conterrà la somma degli elementi del vettore x
    somma = 0;

%   Variabile che conterrà le ciffre soggette a cancellazione degli
%   elementi sommati
    c = 0;

%   Lunghezza del vettore in input
    n = length(x);

%   Ciclo di somma [ vedi wikipedia per dettagli ]
    for i = 1:n
        y = x(i) - c;
        t = somma + y;
        c = (t - somma) - y;
        somma = t;
    end
end