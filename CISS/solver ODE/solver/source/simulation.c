#include "../include/simulation.h"

#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <math.h>

/*
 *  File che conterra' le definizioni di tutte le funzioni e
 *  strutture dati necessarie alla simulazione dell'evoluzione dello stato
 *  del nostro sistema quantistico
 */

//---------------------------------------------------------------------------//

/*
 *  Funzione per l'allocazione dei parametri della simulazione
 */
parametri_simulazione *alloc_parametri_simulazione(double t_fin, 
                                                    uint32_t N,
                                                    uint32_t N2,
                                                    uint32_t n_oss,
                                                    uint64_t dim,
                                                    char *f_av)
{
    parametri_simulazione *param = 
        (parametri_simulazione*)malloc(sizeof(parametri_simulazione));

/*
 *  Verifichiamo l'esito dell'allocazione
 */
    if (param == NULL)
    {
        fprintf(stderr, "ERRORE nell'allocazione dei parametri della "
                        "simulazione.\n"
                        "Memoria insufficiente.\n"
                        "ABORTING...\n");
        return NULL;
    }
/*
 *  Inizializziamo i membri della struttura
 */
    param->t_fin = t_fin;
    param->n_oss = n_oss;
    param->N  = N;
    param->N2 = N2;

/*
 *  Inizializziamo gli stream su file di output
 */
    param->out_Re  = fopen("./oss_RE.bin", "w");
    param->out_Im  = fopen("./oss_IM.bin", "w");
    param->out_rho = fopen("./rho_out.bin", "w");
    param->out_avanzamento = fopen(f_av, "w");
/*
 *  Verifichiamo l'apertura degli stream
 */
    if (param->out_Re == NULL || param->out_Im == NULL 
                              || param->out_rho == NULL
                              || param->out_avanzamento == NULL)
    {
        fprintf(stderr, "ERRORE nell'apertura degli stream su file di "
                        "output per la stampa dei valori di aspettazione "
                        "delle osservabili fornite in input.\n"
                        "Memoria insufficiente.\n"
                        "ABORTING...\n");
        return NULL;
    }
/*
 *  Inizializziamo lo stream su file di output
 */

/*
 *  Allochiamo i membri della struttura
 */
    param->rho = gsl_matrix_complex_calloc(dim, dim);
/*
 *  Verifichiamo l'esito dell'allocazione
 */
    if (param->rho == NULL)
    {
        fprintf(stderr, "ERRORE nell'allocazione dei parametri della "
                        "simulazione.\n"
                        "Memoria insufficiente.\n"
                        "ABORTING...\n");
        free_parametri_simulazione(&param);
        return NULL;
    }
/*
 *  Allochiamo lo spazio per gli indirizzi delle osservabili
 */
    param->oss = 
        (gsl_matrix_complex**)calloc(n_oss, sizeof(gsl_matrix_complex*));
/*
 *  Verifichiamo l'esito dell'allocazione
 */
    if (param->oss == NULL)
    {
        fprintf(stderr, "ERRORE nell'allocazione dei parametri della "
                        "simulazione.\n"
                        "Memoria insufficiente.\n"
                        "ABORTING...\n");
        free_parametri_simulazione(&param);
        return NULL;
    }
/*
 *  Allochiamo lo spazio per tutte le osservabili
 */
    for (uint32_t i = 0; i < n_oss; ++i)
    {
        param->oss[i] = gsl_matrix_complex_calloc(dim, dim);
/*
 *      Verifichiamo l'esito dell'allocazione
 */
        if (param->oss[i] == NULL)
        {
            fprintf(stderr, "ERRORE nell'allocazione dei parametri della "
                            "simulazione.\n"
                            "Memoria insufficiente.\n"
                            "ABORTING...\n");
            free_parametri_simulazione(&param);
            return NULL;
        }
    }

    return param;
}

/*
 *  Funzione per la deallocazione dei parametri della simulazione
 */
void free_parametri_simulazione(parametri_simulazione **param)
{
    if (param == NULL || *param == NULL)
    {
        return;
    }

    if((*param)->rho != NULL)
    {
        gsl_matrix_complex_free((*param)->rho);
        (*param)->rho = NULL;
    }

    for (uint32_t i = 0; i < (*param)->n_oss; ++i)
    {
        if((*param)->oss[i] != NULL)
        {
            gsl_matrix_complex_free((*param)->oss[i]);
            (*param)->oss[i] = NULL;
        }
    }

    if((*param)->out_Re != NULL)
    {
        fclose((*param)->out_Re);
        (*param)->out_Re = NULL;
    }

    if((*param)->out_Im != NULL)
    {
        fclose((*param)->out_Im);
        (*param)->out_Im = NULL;
    }

    if((*param)->out_rho != NULL)
    {
        fclose((*param)->out_rho);
        (*param)->out_rho = NULL;
    }

    if((*param)->out_avanzamento!= NULL)
    {
        fclose((*param)->out_avanzamento);
        (*param)->out_avanzamento = NULL;
    }

    (*param)->t_fin = 0.0;
    (*param)->n_oss = 0;
    (*param)->N = 0;

    free(*param);
    *param = NULL;

    return;
}

//---------------------------------------------------------------------------//

/*
 *  Funzione che calcola le osservabili fornite in input al progamma
 *  all'istante passato e le stampa sull'opportuno file di output
 */
void calcola_osservabili(parametri_simulazione *param,
                            uint32_t dim,
                            double *sol,
                            double t)
{
/*
 *  Variabili d'appoggio
 */
    double pop_0 = 0;
    double pop_1 = 0;
    double oss_0 = 0;
    double oss_1 = 0;
    double oss_2 = 0;
    double oss_3 = 0;
    double oss_4 = 0;
    double oss_p = 0;
    double oss_u = 0;
    double somma = 0;
    double diff  = 0;

/*
 *  Stampiamo su file il valore del tempo attuale
 */
    fwrite(&t, sizeof(double), 1, param->out_Re);
    fwrite(&t, sizeof(double), 1, param->out_Im);

/*
 *  Generiamo una matrix view della soluzione attuale
 */
    gsl_matrix_complex_view rho =
        gsl_matrix_complex_view_array(sol, dim, dim);

/*
 *  Recuperiamo le popolazioni dei primi due stati del sistema
 */
    gsl_complex tmp_0 = gsl_matrix_complex_get(&rho.matrix, 0, 0);
    gsl_complex tmp_1 = gsl_matrix_complex_get(&rho.matrix, 1, 1);

    pop_0 = tmp_0.dat[0];
    pop_1 = tmp_1.dat[0];

/*
 *  Cicliamo su tutte le osservabili da calcolare
 */
    for (uint32_t i = 0; i < param->n_oss; ++i)
    {
/*
 *      Variabile d'appoggio per il calcolo della traccia
 */
        gsl_complex dotu = gsl_complex_rect(0.0, 0.0);
/*
 *      Inizializziamo il valore della traccia
 */
        gsl_complex trace = gsl_complex_rect(0.0, 0.0);
/*
 *      Recuperiamo l'osservabile attuale
 */
        gsl_matrix_complex *oss = param->oss[i];

/*
 *      Calcoliamo la traccia della matrice
 */
        for(uint32_t p = 0; p < dim; ++p)
        {
/*
 *          Generiamo la vector view della p-esima riga di rho
 */
            gsl_vector_complex_view riga = 
                gsl_matrix_complex_row(&rho.matrix, p);
/*
 *          Generiamo la vector view della p-esima colonna dell'osservabile
 */
            gsl_vector_complex_view colonna = 
                gsl_matrix_complex_column(oss, p);
/*
 *          Calcoliamo il risultato opportuno
 */
            gsl_blas_zdotu(&riga.vector, &colonna.vector, &dotu);
/*
 *          Aggiorniamo in modo opportuno il valore della traccia calcolato
 */
            trace.dat[0] = trace.dat[0] + dotu.dat[0];
            trace.dat[1] = trace.dat[1] + dotu.dat[1];
        }

/*
 *      Stampiamo su file il valore d'aspettazione calcolato per
 *      l'osservabile attuale
 */
        fwrite(&(trace.dat[0]), sizeof(double), 1, param->out_Re);
        fwrite(&(trace.dat[1]), sizeof(double), 1, param->out_Im);

/*
 *      Salviamo i valori delle osservabili interessanti
 */
        if (i < 6 || i > param->n_oss - 3)
        {
            switch(i){
                case 0:
                    oss_0 = trace.dat[0];
                case 1:
                    break;
                case 2:
                    oss_1 = trace.dat[0];
                    break;
                case 3:
                    oss_2 = trace.dat[0];
                    break;
                case 4:
                    oss_3 = trace.dat[0];
                    break;
                case 5:
                    oss_4 = trace.dat[0];
                    break;
                default:
                    if (i == param->n_oss - 2){
                        oss_p = trace.dat[0];
                    }
                    else
                    {
                        oss_u = trace.dat[0];
                    }
                        
                    break;
            }
        }
    }

/*
 *  Facciamo un flush del buffer del file delle osservabili per essere sicuri
 *  che tutte le informazioni siano aggiornate
 */
    fflush(param->out_Re);
    fflush(param->out_Im);

/*
 *  Stampiamo la traccia di rho in output
 */
    fprintf(param->out_avanzamento, "%.6e\t", oss_0);

/*
 *  Calcoliamo somma e differenza delle popolazioni dei primi due stati del
 *  sistema
 */
    somma = pop_1 + pop_0;
    diff  = pop_1 - pop_0;
/*
 *  Stampiamo su stdout i valori delle prime due osservabili
 */
    if (somma > 1e-12)
    {
        fprintf(param->out_avanzamento, "%.2e\t%.2e\t%.2e\t", 
                                        diff, somma, diff/somma);
    }
    else
    {
        fprintf(param->out_avanzamento, "%.2e\t%.2e\t%.2e\t", 
                                        diff, somma, 0.0);
    }

/*
 *  Calcoliamo somma e differenza delle prime due osservabili
 */
    somma = oss_1 + oss_2;
    diff  = oss_1 - oss_2;
/*
 *  Stampiamo su stdout i valori delle prime due osservabili
 */
    if (somma > 1e-12)
    {
        fprintf(param->out_avanzamento, "%.2e\t%.2e\t", 
                                        somma, diff/somma);
    }
    else
    {
        fprintf(param->out_avanzamento, "%.2e\t%.2e\t", 
                                        somma, 0.0);
    }

/*
 *  Calcoliamo somma e differenza delle seconde due osservabili
 */
    somma = oss_3 + oss_4;
    diff  = oss_3 - oss_4;
/*
 *  Stampiamo su stdout i valori delle seconde due osservabili
 */
    if (somma > 1e-12)
    {
        fprintf(param->out_avanzamento, "%.2e\t%.2e\t", somma, diff/somma);
    }
    else
    {
        fprintf(param->out_avanzamento, "%.2e\t%.2e\t", somma, 0.0);
    }

/*
 *  Calcoliamo somma e differenza delle ultime due osservabili
 */
    somma = oss_p + oss_u;
    diff  = oss_p - oss_u;
/*
 *  Stampiamo su stdout i valori delle ultime due osservabili
 */
    if (somma > 1e-12)
    {
        fprintf(param->out_avanzamento, "%.2e\t%.2e\n", somma, diff/somma);
    }
    else
    {
        fprintf(param->out_avanzamento, "%.2e\t%.2e\n", somma, 0.0);
    }

/*
 *  Facciamo un flush dello stream su file per essere sicuri che le
 *  informazioni siano sempre aggiornate
 */
    fflush(param->out_avanzamento);

    return;
}

/*
 *  Funzione che stampa sul file di output opportuno la matrice di densità del
 *  sistema
 */
void stampa_rho(parametri_simulazione *param, 
                uint32_t dim, 
                double *sol)
{
/*
 *  Creazione di una matrix view della nostra rho
 */
    gsl_matrix_complex_view rho =
        gsl_matrix_complex_view_array(sol, dim, dim);

/*
 *  DBG richiesto da alessandro.
 *
 *  Stampiamo in output tutta la rho ad ogni istante temporale
 */
    gsl_matrix_complex_fwrite(param->out_rho, &rho.matrix);

/*
 *  Creiamo una vector view della diagonale di rho
 */
    // gsl_vector_complex_view diag = 
        // gsl_matrix_complex_diagonal(&rho.matrix);

/*
 *  Stampiamo su stdout la diagonale di rho
 */
    // gsl_vector_complex_fwrite(param->out_rho, &diag.vector);

/*
 *  Facciamo un flush del buffer per essere sicuri che le informazioni
 *  sul file siano aggiornate
 */
    fflush(param->out_rho);

    return;
}

/*
 *  Funzione che stampa polarizzazione e carica traferita sul donore al
 *  termine della simulazione 
 */
void stampa_risultati_finali(parametri_simulazione *param, 
                                uint32_t dim,
                                double *sol,
                                FILE *out)
{
/*
 *  Variabili d'appoggio
 */
    double oss_0 = 0;
    double oss_1 = 0;

/*
 *  Creazione di una matrix view della nostra rho
 */
    gsl_matrix_complex_view rho =
        gsl_matrix_complex_view_array(sol, dim, dim);

/*
 *  Cicliamo sulle osservabili da calcolare (le ultime due)
 */
    for (uint32_t i = param->n_oss-2; i < param->n_oss; ++i)
    {
/*
 *      Variabile d'appoggio per il calcolo della traccia
 */
        gsl_complex dotu = gsl_complex_rect(0.0, 0.0);
/*
 *      Inizializziamo il valore della traccia
 */
        gsl_complex trace = gsl_complex_rect(0.0, 0.0);
/*
 *      Recuperiamo l'osservabile attuale
 */
        gsl_matrix_complex *oss = param->oss[i];

/*
 *      Calcoliamo la traccia della matrice
 */
        for(uint32_t p = 0; p < dim; ++p)
        {
/*
 *          Generiamo la vector view della p-esima riga di rho
 */
            gsl_vector_complex_view riga = 
                gsl_matrix_complex_row(&rho.matrix, p);
/*
 *          Generiamo la vector view della p-esima colonna dell'osservabile
 */
            gsl_vector_complex_view colonna = 
                gsl_matrix_complex_column(oss, p);
/*
 *          Calcoliamo il risultato opportuno
 */
            gsl_blas_zdotu(&riga.vector, &colonna.vector, &dotu);
/*
 *          Aggiorniamo in modo opportuno il valore della traccia calcolato
 */
            trace.dat[0] = trace.dat[0] + dotu.dat[0];
            trace.dat[1] = trace.dat[1] + dotu.dat[1];
        }

/*
 *      Stampiamo su file il valore d'aspettazione calcolato per
 *      l'osservabile attuale
 */
        fwrite(&(trace.dat[0]), sizeof(double), 1, param->out_Re);
        fwrite(&(trace.dat[1]), sizeof(double), 1, param->out_Im);

/*
 *      Salviamo i valori delle osservabili interessanti
 */
        if (i == param->n_oss-2)
        {
            oss_0 = trace.dat[0];
        }
        else
        {
            oss_1 = trace.dat[0];
        }
    }

/*
 *  Calcoliamo somma e differenza delle osservabili
 */
    double somma = oss_0 + oss_1;
    double diff  = oss_0 - oss_1;
/*
 *  Stampiamo i valori ottenuti
 */
    if (somma > 1e-12)
    {
        fprintf(out, "%.8e\t%.8e\t", somma, diff/somma);
    }
    else
    {
        fprintf(out, "%.8e\t%.8e\t", somma, 0.0);
    }

    return;
}

/*
 *  Funzione che esegue il controllo sulla positività degli autovalori di rho
 */
extern
uint8_t controllo_aut_rho(double *vec, uint64_t dim, double tol)
{
/*
 *  Valore di ritorno della funzione
 */
    uint8_t ris = 1;

/*
 *  Creazione di una matrix view della nostra rho
 */
    gsl_matrix_complex_view rho =
        gsl_matrix_complex_view_array(vec, dim, dim);

/*
 *  Allocazione del workspace per la diagonalizzazione di rho
 */
    gsl_eigen_herm_workspace *wrk = gsl_eigen_herm_alloc(dim);

/*
 *  Verifichiamo l'esito dell'allocazione
 */
    if (wrk == NULL)
    {
        fflush(stderr);
        return 2;
    }

/*
 *  Allochiamo il vettore che conterrà gli autovalori di rho
 */
    gsl_vector *autovalori = gsl_vector_alloc(dim);

/*
 *  Verifichiamo l'esito dell'allocazione
 */
    if (autovalori == NULL)
    {
        gsl_eigen_herm_free(wrk);
        fflush(stderr);
        return 2;
    }

/*
 *  Eseguiamo la diagonalizzazione di rho.
 *
 *  ATTENZIONE:: il processo è distruttivo e rho verrà modificata.
 */
    gsl_eigen_herm(&rho.matrix, autovalori, wrk);

/*
 *  Recuperiamo il minimo autovalore di rho
 */
    double min = gsl_vector_min(autovalori);

/*
 *  Eseguiamo il test su questo autovalore
 */
    if (min + tol < 0)
    {
        ris = 0;
    }

/*
 *  Deallocazione del workspace per la diagonalizzazione di rho
 */
    gsl_eigen_herm_free(wrk);
    gsl_vector_free(autovalori);

    return ris;
}

/*
 *  Funzione che approssima l'evoluzione temporale del sistema considerato
 */
void simula(parametri_simulazione *param, spazio_lavoro_forzante *slf)
{
/*
 *  CHECK su formato input parametri
 */
    // FILE *f_check = fopen("./CHECK_input.bin", "w");
    // gsl_matrix_complex_fwrite(f_check, param->oss[0]);
    // fclose(f_check);

/*
 *  Output avanzamento simulazione
 */
    fprintf(param->out_avanzamento, "Inizio simulazione\n");
/*
 *  Tempo finale della simulazione
 */
    double T = param->t_fin;

/*
 *  Passo temporale a cui stampare le varie osservabili
 */
    double dT = T / (double)param->N;

/*
 *  Contatore necessario per calcolare ogni quanto stampare la rho del sistema
 */
    uint32_t count_rho = 0;
/*
 *  Calcoliamo la dimensione del sistema
 */
    size_t dim  = slf->dim;
    size_t size = dim * dim;
/* 
 *  Spezziamo le condizioni iniziali complesse in parte reale
 *  e parte immaginaria.
 * 
 *  Allochiamo dinamicamente questo vettore per evitare problemi se la
 *  dimensione diventa troppo grossa
 */
    double *ci = (double *)calloc(2 * size, sizeof(double));

/*
 *  Verifichiamo l'esito dell'allocazione
 */
    if (ci == NULL)
    {
        fprintf(stderr, "Errore nell'allocazione della memoria per il vettore "
                        "contennete lo stato del sistema.\n"
                        "ABORTING...\n\n");
        return;
    }

/*
 *  IMPORTANTE: Siccome poi utilizzeremo delle matrix view di questo array
 *              e' necesario generare le condizioni iniziali nell'ordine
 *              corretto. Ovvero, scorre la matrice rho prima per righe
 *              poi per colonne. In altre parole, si fissa la riga e si
 *              scorrono le colonne.
 *              Inoltre dobbiamo salvare i numeri complessi come coppia
 *              consecutia di numeri reali.
 */
    for (size_t i = 0; i < dim; ++i)
    {
        size_t base = i * 2 * dim;
        for (size_t j = 0, k = 0; j < dim; ++j, k += 2)
        {
            gsl_complex tmp = gsl_matrix_complex_get(param->rho, i, j);
            ci[base + k]     = tmp.dat[0];
            ci[base + k + 1] = tmp.dat[1];
        }
    }

/*
 *  Output avanzamento simulazione
 */
    fprintf(param->out_avanzamento, "%.2f\t", 0.0);

/* 
 *  Stampiamo la diagonale di rho all'istante inizale della simulazione
 */
    stampa_rho(param, dim, ci);

/*
 *  Calcoliamo e stampiamo gli osservabili all'istante iniziale della
 *  simulazione
 */
    calcola_osservabili(param, dim, ci, 0.0);

/*
 *  Inizializziamo l'ambiente richiesto dalle interfacce GSL per la
 *  risoluzione del nostro sistema di ODE
 */
    gsl_odeiv2_system s = {forzante, jacobiano, 2*size, slf};

/*
 *  Passo iniziale di integrazione
 */
    // double passo = dT < 1 ? 1e-14 : dT * 1e-14;
    double passo = 5e-4;

/*
 *  Settiamo l'errore target del solver di ode selezionato
 */
    double prec = 1e-10;

/*
 *  Ciclo per l'approssimazione della soluzione nei vari intervalli
 *  richiesti
 */
    for (size_t i = 0; i < param->N; ++i)
    {
/*
 *      Tempo iniziale step di integrazione
 */
        double t_i = i * dT;
/*
 *      Tempo finale step di integrazione
 */
        double t_f = (i + 1) * dT;
/*
 *      Allochiamo il solver GLS 
 *
 *      In ordine:
 *          sistema
 *          metodo numerico, in questo caso il metodo predictor
 *          corrector (Adams-Bashforth) (Adams-Moulton)
 *          passo
 *          errore_assoluto (possibile modificare per aggiustare la precisione)
 *          errore_relativo (possibile modificare per aggiustare la precisione)
 *          peso_errore_funzione
 *          peso_errore_derivata
 */
        gsl_odeiv2_driver *solver_ode =
            gsl_odeiv2_driver_alloc_standard_new(&s,
                gsl_odeiv2_step_rkck, passo, prec, prec, 1.0, 1.0);
/*
 *      Applichiamo il solver GLS per l'intervallo temporale (t_i,t_f)
 */
#ifndef DEBUG   
        gsl_odeiv2_driver_apply(solver_ode, &t_i, t_f, ci);
#else
        int ris = gsl_odeiv2_driver_apply(solver_ode, &t_i, t_f, ci);
        if (ris != GSL_SUCCESS)
        {
            if (ris == GSL_EBADFUNC)
            {
                fprintf(stderr, "Possibile richiamo dello jacobiano da "
                                "parte del solver\n\n");
            }

            fprintf(stderr, "ERRORE nell'evoluzione temporale del "
                            "sistema\n");
            fprintf(stderr, "CODICE ERRORE: %d\n\n", ris);
            fprintf(stderr, "ABORTING...\n");
            return;
        }
#endif

/*
 *      Incrementiamo il contatore che tiene traccia di ogni quanto abbiamo
 *      la necessità di stampare la rho del sistema
 */
            ++count_rho;

/*
 *      Controlliamo se dobbiamo stampare o meno la rho
 */
        if (!(count_rho < param->N2))
        {
/*
 *          Stampiamo su stdout la diagonale della matrice di densità del
 *          sistema
 */
            stampa_rho(param, dim, ci);
/*
 *          Reinizializziamo il contatore
 */
            count_rho = 0;
        }

/*
 *      Output avanzamento simulazione
 */
        fprintf(param->out_avanzamento, "%.2f\t", t_f);

/*
 *      Calcoliamo e stampiamo gli osservabili all'istante attuale della
 *      simulazione
 */
        calcola_osservabili(param, dim, ci, t_f);

/*
 *      Deallochiamo il solverl GSL per evitare un utilizzo eccessivo della RAM
 */
        gsl_odeiv2_driver_free(solver_ode);
    }

/*
 *  Stampa polarizazzione e carica trasferita finali
 */
    FILE *f_ris_fin = fopen("./final_results.txt","w");
    stampa_risultati_finali(param, dim, ci, f_ris_fin);
    fclose(f_ris_fin);

/*
 *  Creazione di una matrix view della nostra rho
 */
    gsl_matrix_complex_view rho =
        gsl_matrix_complex_view_array(ci, dim, dim);

/*
 *  Stampa matrice di densità finale
 */
    FILE *f_rho_fin = fopen("./rho_fin.bin", "w");
    gsl_matrix_complex_fwrite(f_rho_fin, &rho.matrix);
    fclose(f_rho_fin);

/*
 *  Aggiornamento stato avanzamento simulazione
 */
    fprintf(param->out_avanzamento, "TEST AUTOVALORI:: Inizio.\n");
    fflush(param->out_avanzamento);

/*
 *  Eseguiamo il controllo sulla positività degli autovalori di rho
 */
    uint8_t ris = controllo_aut_rho(ci, dim, 10 * prec * param->N);

/*
 *  Verifichiamo l'esito del controllo della positività degli autovalori di rho
 */
    if (ris == 1)
    {
        fprintf(param->out_avanzamento, "TEST AUTOVALORI:: OK\n");
        fflush(param->out_avanzamento);
    } 
    else if (ris == 0)
    {
        fprintf(param->out_avanzamento, "TEST AUTOVALORI:: X\n"
                                        "ATTENZIONE:: autovalori di rho "
                                        "negativi oltre la precisione attesa"
                                        "del risolutore numerico del problema "
                                        "differenziale.\n");
        fflush(param->out_avanzamento);
    }
    else
    {
        fprintf(param->out_avanzamento, "TEST AUTOVALORI:: X\n"
                                        "ATTENZIONE:: allocazione memoria per "
                                        "il calcolo degli autovalori di rho "
                                        "non riuscita.\n");
        fflush(param->out_avanzamento);
    }

/*
 *  Deallochiamo la memoria riservta al vettore contenente lo stato del sistema
 */
    free(ci);

    return;
}