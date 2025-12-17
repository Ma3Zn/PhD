#ifndef EVOLUTION_H
#define EVOLUTION_H

#include <gsl/gsl_matrix.h>
#include <stdint.h>

/*
 *  Header file contenente le dichiarazioni di tutte le funzioni e
 *  strutture dati necessarie alle interfacce GSL per poter eseguire 
 *  l'approssimazione numerica del sistema di equazioni differenziali
 *  che governano lo stato del sistema quantistico
 */

//---------------------------------------------------------------------------//

/*
 *  Struttura dati contenente tutte le variabili d'appoggio necessarie per il
 *  calcolo del forzante
 */
typedef struct
{
/*
 *  Dimensione dello spazio di Hilbert considerato.
 */
    uint32_t dim;
/*
 *  Contiene la matrice di densità del sistema durante l'evoluzione temporale
 */
    gsl_matrix_complex *rho;
/*
 *  Contiene l'Hamiltoniana tempo indipendente del sistema
 */
    gsl_matrix_complex *H;
/*
 *  Numero di operatori d'errore
 */
    uint32_t n_err;
/*
 *  Array di puntatori agli operatori d'errore
 */
    gsl_matrix_complex **operatore_errore;
/*
 *  Array di puntatori agli operatori d'errore fittizi della nuova forma 
 *  dell'equazione master
 */
    gsl_matrix_complex **operatore_errore_fittizio;
/*
 *  Array contenente i ratei di errore
 */
    double *rateo;
/*
 *  Variabile d'appoggio per il calcolo del termine incoerente della master
 *  equation
 */
    gsl_matrix_complex *tmp;
} spazio_lavoro_forzante;

//---------------------------------------------------------------------------//

/*
 *  Funzione per l'allocazione dello spazio di lavoro per il calcolo del
 *  termine forzante
 */
extern
spazio_lavoro_forzante *alloc_spazio_lavoro_forzante(uint32_t dim, 
                                                        uint32_t n_err);

/*
 *  Funzione per la deallocazione dello spazio di lavoro per il calcolo del
 *  termine forzante
 */
extern
void free_spazio_lavoro_forzante(spazio_lavoro_forzante **slf);

//---------------------------------------------------------------------------//

/* 
 *  Dichiarazione delle due funzioni che permettono di valutare il sistema
 *  di ODE tramite le procedure fornite da GSL
 */
extern int forzante(double t, const double y[], double f[], void *params);
extern int jacobiano(double t, const double y[], double *dfdy, double dfdt[],
                        void *params);

/*
 *  Funzione che aggiunge l'errore secondo la Limblad Master Equation
 */
extern
uint8_t aggiungi_errore(gsl_matrix_complex *rho, spazio_lavoro_forzante *slf);
#endif