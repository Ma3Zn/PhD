#ifndef SIMULATION_H
#define SIMULATION_H

#include "evolution.h"

#include <stdbool.h>
#include <stdio.h>

/*
 *  Header file che conterra' le dichiarazione di tutte le funzioni e
 *  strutture dati necessarie alla simulazione dell'evoluzione dello stato
 *  del nostro sistema quantistico
 */

//---------------------------------------------------------------------------//

/*
 *  Struttura dati contenente tutti i parametri necessari alla simulazione
 */
typedef struct
{
/*
 *  Tempo finale della simulazione
 */
    double t_fin;
/*
 *  Numero di intervalli in cui decomporre la mesh temporale dell'evoluzione
 */
    uint32_t N;
/*
 *  Numero di istanti temporali in cui stampare la rho completa del sistema
 */
    uint32_t N2;
/*
 *  Condizione iniziale del sistema
 */
    gsl_matrix_complex *rho;
/*
 *  Numero di osservabili
 */
    uint32_t n_oss;
/*
 *  Osservabili del sistema
 */
    gsl_matrix_complex **oss;
/*
 *  Stream file di output (parte reale dei valori di aspettazione)
 */
    FILE *out_Re;
/*
 *  Stream file di output (parte immaginario dei valori di aspettazione)
 */
    FILE *out_Im;
/*
 *  Stream file di output per la matrice di densità del sistema
 */
    FILE *out_rho;
/*
 *  Stream file di output in cui stampare l'avanzamento della simulazione
 */
    FILE *out_avanzamento;
}parametri_simulazione;

//---------------------------------------------------------------------------//

/*
 *  Funzione per l'allocazione dei parametri della simulazione
 */
extern
parametri_simulazione *alloc_parametri_simulazione(double t_fin, 
                                                    uint32_t N,
                                                    uint32_t N2,
                                                    uint32_t n_oss,
                                                    uint64_t dim,
                                                    char *f_av);

/*
 *  Funzione per la deallocazione dei parametri della simulazione
 */
extern
void free_parametri_simulazione(parametri_simulazione **param);

//---------------------------------------------------------------------------//

/*
 *  Funzione che calcola le osservabili fornite in input al progamma
 *  all'istante passato e le stampa sull'opportuno file di output
 */
extern
void calcola_osservabili(parametri_simulazione *param,
                            uint32_t dim,
                            double *sol,
                            double t);

/*
 *  Funzione che stampa sul file di output opportuno la matrice di densità del
 *  sistema
 */
extern
void stampa_rho(parametri_simulazione *param, 
                uint32_t dim, 
                double *sol);

/*
 *  Funzione che stampa polarizzazione e carica traferita sul donore al
 *  termine della simulazione 
 */
extern
void stampa_risultati_finali(parametri_simulazione *param, 
                                uint32_t dim,
                                double *sol,
                                FILE *out);

/*
 *  Funzione che esegue il controllo sulla positività degli autovalori di rho
 */
extern
uint8_t controllo_aut_rho(double *vec, uint64_t dim, double tol);

/*
 *  Funzione che approssima l'evoluzione temporale del sistema considerato
 */
extern
void simula(parametri_simulazione *param, spazio_lavoro_forzante *slf);

#endif