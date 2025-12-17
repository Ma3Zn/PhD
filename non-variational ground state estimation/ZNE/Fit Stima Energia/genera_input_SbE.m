function SbE = genera_input_SbE(dir, tstep, ordine, data_range, prec_est, N_sample)
%   Funzione che recupera tutti dati necessari per eseguire la SbE con ZNE 
%   all'ordine opportuno

%%  Recuperiamo i valori di aspettazione per la ZNE

%   Individuiamo il percorso opportuno ai dati
    dir             = strcat(dir, "/", num2str(tstep));
    exp_val         = load(strcat(dir, "/exp_values.txt"));
    folding_factors = load(strcat(dir, "/folding_factors.txt"));

%   Calcoliamo quante colonne recuperare dai file per eseguire una SbE di
%   un dato ordine e filtriamo i fati per gli opportuni folding factors.
%   Ne facciamo poi la trasposizione per riportare i dati nel formato
%   richiesto dalle interfacce ai metodi di fit dei dati
    pow_max         = 2 * ordine + 1;

    exp_val         = exp_val(data_range, 1:pow_max).';
    folding_factors = folding_factors(data_range).';

%   Aggiungiamo una varianza ai dati
    ZNE.exp_val = aggiungi_varianza(exp_val, prec_est, N_sample);
    ZNE.ff      = aumenta_ff(folding_factors, N_sample);
   
%   Salviamo i dati in una struct opportuna
    SbE.folding_factors = folding_factors;
    SbE.data_range      = data_range;
    SbE.prec_est        = prec_est;
    SbE.ZNE             = ZNE;
end