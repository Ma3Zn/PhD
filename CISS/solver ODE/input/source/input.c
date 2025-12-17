#include "../include/input.h"

/*
 *  Header file contenente le definizioni di tutte le funzioni per
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
                    char *file_oss)
{
/*
 *  Apriamo gli stream di input
 */
    FILE *f_rho = fopen(file_rho_in, "r");
    FILE *f_H   = fopen(file_H, "r");
    FILE *f_err = fopen(file_err, "r");
    FILE *f_rat = fopen(file_ratei, "r");
    FILE *f_oss = fopen(file_oss, "r");

/*
 *  Verifichiamo l'apertura dei file
 */
    if (f_rho == NULL)
    {
        fprintf(stderr, "Errore nell'apertura del file di input di rho.\n"
                        "ABORTING...\n\n");
        return 1;
    }

    if (f_H == NULL)
    {
        fprintf(stderr, "Errore nell'apertura del file di input di H.\n"
                        "ABORTING...\n\n");
        fclose(f_rho);
        return 1;
    }

    if (f_err == NULL)
    {
        fprintf(stderr, "Errore nell'apertura del file di input degli errori."
                        "\n"
                        "ABORTING...\n\n");
        fclose(f_rho);
        fclose(f_H);
        return 1;
    }

    if (f_rat == NULL)
    {
        fprintf(stderr, "Errore nell'apertura del file di input delle "
                        "osservabili.\n"
                        "ABORTING...\n\n");
        fclose(f_rho);
        fclose(f_H);
        fclose(f_err);
        return 1;
    }

    if (f_oss == NULL)
    {
        fprintf(stderr, "Errore nell'apertura del file di input delle "
                        "osservabili.\n"
                        "ABORTING...\n\n");
        fclose(f_rho);
        fclose(f_H);
        fclose(f_err);
        fclose(f_rat);
        return 1;
    }

/*
 *  Variabile che conterrà l'esito della lettura dei vari operatori da file
 */
    size_t ris;
/*
 *  Leggiamo rho_in
 */
    ris = gsl_matrix_complex_fread(f_rho, param->rho);
/*
 *  Verifichiamo l'esito della lettura
 */
    if (ris == GSL_EFAILED)
    {
        fprintf(stderr, "Errore nella lettura di rho.\n"
                        "ABORTING...\n\n");
        fclose(f_rho);
        fclose(f_H);
        fclose(f_err);
        fclose(f_rat);
        fclose(f_oss);

        return 1;
    }
/*
 *  Leggiamo H
 */
    ris = gsl_matrix_complex_fread(f_H, slf->H);
/*
 *  Verifichiamo l'esito della lettura
 */
    if (ris == GSL_EFAILED)
    {
        fprintf(stderr, "Errore nella lettura di H.\n"
                        "ABORTING...\n\n");
        fclose(f_rho);
        fclose(f_H);
        fclose(f_err);
        fclose(f_rat);
        fclose(f_oss);

        return 1;
    }
/*
 *  Leggiamo tutti gli operatori d'errore
 */
    for (uint32_t i = 0; i < slf->n_err; ++i)
    {
        ris = gsl_matrix_complex_fread(f_err, slf->operatore_errore[i]);
/*
 *      Verifichiamo l'esito della lettura
 */
        if (ris == GSL_EFAILED)
        {
            fprintf(stderr, "Errore nella lettura degli errori.\n"
                            "ABORTING...\n\n");
            fclose(f_rho);
            fclose(f_H);
            fclose(f_err);
            fclose(f_rat);
            fclose(f_oss);

            return 1;
        }
    }
/*
 *  Leggiamo tutti gli operatori d'errore fittizi che sono riportati nel file
 *  degli operatori d'errore dopo gli operatiri d'errore standard
 */
    for (uint32_t i = 0; i < slf->n_err; ++i)
    {
        ris = gsl_matrix_complex_fread(f_err, slf->operatore_errore_fittizio[i]);
/*
 *      Verifichiamo l'esito della lettura
 */
        if (ris == GSL_EFAILED)
        {
            fprintf(stderr, "Errore nella lettura degli errori.\n"
                            "ABORTING...\n\n");
            fclose(f_rho);
            fclose(f_H);
            fclose(f_err);
            fclose(f_rat);
            fclose(f_oss);

            return 1;
        }
    }
/*
 *  leggiamo i ratei opportuno
 */
    ris = fread(slf->rateo, sizeof(double), slf->n_err, f_rat);
/*
 *  Verifichiamo l'esito della lettura
 */
    if (ris != slf->n_err)
    {
        fprintf(stderr, "Errore nella lettura dei ratei d'errore.\n"
                            "ABORTING...\n\n");
        fclose(f_rho);
        fclose(f_H);
        fclose(f_err);
        fclose(f_rat);
        fclose(f_oss);

        return 1;
    }
/*
 *  Leggiamo tutte le osservabili
 */
    for (uint32_t i = 0; i < param->n_oss; ++i)
    {
        ris = gsl_matrix_complex_fread(f_oss, param->oss[i]);
/*
 *      Verifichiamo l'esito della lettura
 */
        if (ris == GSL_EFAILED)
        {
            fprintf(stderr, "Errore nella lettura delle osservabili.\n"
                            "ABORTING...\n\n");
            fclose(f_rho);
            fclose(f_H);
            fclose(f_err);
            fclose(f_rat);
            fclose(f_oss);

            return 1;
        }
    }

/*
 *  Chiudiamo tutti gli stream di input
 */
    fclose(f_rho);
    fclose(f_H);
    fclose(f_err);
    fclose(f_rat);
    fclose(f_oss);

    return 0;
}