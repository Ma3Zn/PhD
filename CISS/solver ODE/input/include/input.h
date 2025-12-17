#ifndef INPUT_H
#define INPUT_H

#include <stdio.h>

#include "../../solver/include/evolution.h"
#include "../../solver/include/simulation.h"

/*
 *  Header file contenente le dichiarazioni di tutte le funzioni per
 *  eseguire la lettra degli operatori passati in input al programma
 */

//---------------------------------------------------------------------------//

/*
 *  Funzione che legge i vari operatori forniti in input
 */
extern
uint8_t leggi_input(parametri_simulazione *param, 
                    spazio_lavoro_forzante *slf,
                    char *file_rho_in, 
                    char *file_H,      
                    char *file_err,
                    char *file_ratei,    
                    char *file_oss);

#endif