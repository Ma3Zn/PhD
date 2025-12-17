#include "../include/evolution.h"

#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_blas.h>

/* 
 *  File che contiene la definizione di tutte le funzioni necessarie
 *  alle prcedure  presenti in GSL per la risoluzione di un sistema di ODE
 */

//---------------------------------------------------------------------------//

/*
 *  Funzione per l'alloczione dello spazio di lavoro per il calcolo del
 *  termine forzante
 */
spazio_lavoro_forzante *alloc_spazio_lavoro_forzante(uint32_t dim, 
                                                            uint32_t n_err)
{
    spazio_lavoro_forzante *slf = 
        (spazio_lavoro_forzante*)malloc(sizeof(spazio_lavoro_forzante));

/*
 *  Verifichiamo l'esito dell'allocazione
 */
    if (slf == NULL)
    {
        fprintf(stderr, "ERRORE nell'allocazione dello spazio di lavoro "
                        "per il forzante.\n"
                        "Memoria insufficiente.\n"
                        "ABORTING...\n");
        return NULL;
    }

/*
 *  Inizializziamo i membri della struttura
 */
    slf->n_err = n_err;
    slf->dim = dim;

/*
 *  Allochiamo i membri della struttura
 */
    slf->rho = gsl_matrix_complex_calloc(dim, dim);
    slf->H   = gsl_matrix_complex_calloc(dim, dim);

    slf->operatore_errore =
        (gsl_matrix_complex**)calloc(n_err, sizeof(gsl_matrix_complex*));
    slf->operatore_errore_fittizio = 
        (gsl_matrix_complex**)calloc(n_err, sizeof(gsl_matrix_complex*));
    slf->rateo = (double *)calloc(n_err, sizeof(double));

/*
 *  Verifichiamo gli esisti delle allocazioni
 */
    if (   slf->rho == NULL || slf->H == NULL || slf->operatore_errore == NULL
        || slf->operatore_errore == NULL || slf->rateo == NULL)
    {
        fprintf(stderr, "ERRORE nell'allocazione dello spazio di lavoro "
                        "per il forzante.\n"
                        "Memoria insufficiente.\n"
                        "ABORTING...\n");
        free_spazio_lavoro_forzante(&slf);
        return NULL;
    }

/*
 *  Allocazione memoria per operatori d'errore
 */
    for (uint32_t i = 0; i < n_err; ++i)
    {
        slf->operatore_errore[i] = gsl_matrix_complex_calloc(dim, dim);
        slf->operatore_errore_fittizio[i] = 
            gsl_matrix_complex_calloc(dim, dim);
/*
 *      Verifichiamo l'esito dell'allocazione
 */
        if (   slf->operatore_errore[i] == NULL 
            || slf->operatore_errore_fittizio[i] == NULL)
        {
            fprintf(stderr, "ERRORE nell'allocazione dello spazio di lavoro "
                            "per il forzante.\n"
                            "Memoria insufficiente.\n"
                            "ABORTING...\n");
            free_spazio_lavoro_forzante(&slf);
            return NULL;
        }
    }

/*
 *  Allocazione memoria variabile d'appoggio per il calcolo del termine
 *  incoerente della master equation
 */
    slf->tmp = gsl_matrix_complex_calloc(dim, dim);
/*
 *  Verifichiamo l'esito dell'allocazione
 */
    if (   slf->tmp == NULL)
    {
        fprintf(stderr, "ERRORE nell'allocazione dello spazio di lavoro "
                        "per il forzante.\n"
                        "Memoria insufficiente.\n"
                        "ABORTING...\n");
        free_spazio_lavoro_forzante(&slf);
        return NULL;
    }

    return slf;    
}

/*
 *  Funzione per la dealloczione dello spazio di lavoro per il calcolo del
 *  termine forzante
 */
void free_spazio_lavoro_forzante(spazio_lavoro_forzante **slf)
{
    if (slf == NULL || *slf == NULL)
    {
        return;
    }

    if ((*slf)->rho != NULL)
    {
        gsl_matrix_complex_free((*slf)->rho);
    }

    if ((*slf)->H != NULL)
    {
        gsl_matrix_complex_free((*slf)->H);
    }

    if ((*slf)->rateo != NULL)
    {
        free((*slf)->rateo);
    }

    if ((*slf)->operatore_errore != NULL)
    {
        for (uint32_t i = 0; i < (*slf)->n_err; ++i)
        {
            gsl_matrix_complex_free((*slf)->operatore_errore[i]);
            gsl_matrix_complex_free((*slf)->operatore_errore_fittizio[i]);
        }
    }

    free((*slf)->operatore_errore);
    free((*slf)->operatore_errore_fittizio);

    (*slf)->n_err = 0;
    (*slf)->dim = 0;

    free(*slf);
    *slf = NULL;

    return;
}

//---------------------------------------------------------------------------//

/* 
 *  Funzione per il calcolo del forzante della master equation
 */
int forzante(double t, const double y[], double f[], void *params)
{
    spazio_lavoro_forzante *slf = (spazio_lavoro_forzante *)params;
    
    size_t dim  = slf->dim;
    size_t size = dim * dim;

/*
 *  Variabili d'appoggio per il calcolo del forzante
 */
    gsl_matrix_complex *rho = slf->rho;
    gsl_matrix_complex *tmp = slf->tmp;

/*
 *  Hamiltonina tempo indipendente del sistema
 */
    gsl_matrix_complex *H = slf->H;

/*
 *  Inizializzazione fattori comuni
 */
    gsl_complex imaginary_unit  = gsl_complex_rect( 0.0, 1.0);
    gsl_complex minus_one       = gsl_complex_rect(-1.0, 0.0);
    gsl_complex zero            = gsl_complex_rect( 0.0, 0.0);
    gsl_complex one             = gsl_complex_rect( 1.0, 0.0);

/*
 *  Vogliamo creare una matrix view del dato iniziale, ma siccome y NON deve
 *  essere modificato, utiliziamo l'array di output f. Dunque come prima cosa
 *  inizializziamo f = y, ricordando che y contiene 2n entrate contententi
 *  coppie di numeri immaginari (Re{z},Im{z})
 * 
 *  ATTEZIONE:: PARALLELIZZABILE!!!!  (Da fare assolutamente)
 */
    for (size_t i = 0; i < 2 * size; ++i)
    {
/*
 *      Copia delle parti reali della soluzione
 */
        f[i] = y[i];
    }

/*
 *  Creazione di una matrix view di f 
 */
    gsl_matrix_complex_view forzante =
        gsl_matrix_complex_view_array(f, dim, dim);

/*
 *  Inizializzazione rho (condizione iniziale del problema)
 */
    gsl_matrix_complex_memcpy(rho, &forzante.matrix);

/*
 *  Calcoliamo il commutatore -i*[H,rho] = i*[rho,H]
 *
 *  termine: i * H * rho
 */
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, imaginary_unit, 
                    H, rho, zero, &forzante.matrix);
/*
 *  termine: i * rho * H - commutatore
 */
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, imaginary_unit,
                    rho, H, minus_one, &forzante.matrix);

/*
 *  Aggiungiamo l'errore
 */
for (size_t i = 0; i < slf->n_err; ++i)
    {
        gsl_matrix_complex *X_i = slf->operatore_errore[i];
        gsl_matrix_complex *A_i = slf->operatore_errore_fittizio[i];
        gsl_complex rateo = gsl_complex_rect(slf->rateo[i], 0.0);

/*
 *      Calcoliamo il termine della master equation
 *                      rateo[i]*X*rho*A
 */
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans,
                        rateo, X_i, rho,
                        zero, tmp);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans,
                        one, tmp, A_i,
                        one, &forzante.matrix);

/*
 *      Calcoliamo il termine della master equation
 *                      -rateo[i]*rho*A*X
 */
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans,
                        rateo, rho, A_i,
                        zero, tmp);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans,
                        minus_one, tmp, X_i,
                        one, &forzante.matrix);

/*
 *      Aggiungiamo ora i termini dagati della master equation
 */


/*
 *      Calcoliamo il termine della master equation
 *                      rateo[i]*A'*rho*X'
 */
        gsl_blas_zgemm(CblasConjTrans, CblasNoTrans,
                        rateo, A_i, rho,
                        zero, tmp);
        gsl_blas_zgemm(CblasNoTrans, CblasConjTrans,
                        one, tmp, X_i,
                        one, &forzante.matrix);

/*
 *      Calcoliamo il termine della master equation
 *                      -rateo[i]*X'*A'*rho
 */
        gsl_blas_zgemm(CblasConjTrans, CblasConjTrans,
                        rateo, X_i, A_i,
                        zero, tmp);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans,
                        minus_one, tmp, rho,
                        one, &forzante.matrix);
    }

    return GSL_SUCCESS;
}

/* 
 * Lo Jacobiano è richiesto dalle interfacce GSL ma non tutti i solver lo
 * utilizzano. In particolare il solver da noi utilizzato non ne fa uso.
 */
int jacobiano(double t, const double y[], double *dfdy, double dfdt[],
                void *params)
{
    return GSL_EBADFUNC;
}