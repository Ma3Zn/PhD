function sim = rigenera_campioni_noisy(sim, data, ordine_SbE, tstep, N_sample)
%   Funzione che rigenera il campione opportuno dei dati da utilizzare per
%   la ZNE

%   Calcoliamo quante righe recuperare dai file per eseguire una SbE di un
%   dato ordine
    pow_max = 2 * ordine_SbE + 1;

%   Recuperiamo il file opportuno da leggere
    file = strcat(data, "/", num2str(tstep), "/exp_values.txt");

%   Carichiamo il file
    exp_val = load(file);

%   Filtriamo i valori di exp_val
    exp_val = exp_val(sim.data_range,1:pow_max).';

%   Eseguiamo il campionamento dei dati
%   ATTENZIONE:: la deviazione standard è la radice della varianza
    sim.ZNE.exp_val = aggiungi_varianza(exp_val, sim.prec_est^(1/2), N_sample);
    sim.ZNE.ff      = aumenta_ff(sim.folding_factors, N_sample);
end