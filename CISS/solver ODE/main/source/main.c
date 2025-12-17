#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#include "../../input/include/input.h"

/*
 *  Entry point
 */
int main(int argc, char **argv)
{
    --argc;
    ++argv;

/*
 *  Verifichiamo il corretto numero di input al programma
 */
    if (argc < 12)
    {
        fprintf(stderr, "Numero di parametri di inpur forniti errato.\n"
                        "In ordine:\n\t1. dimensione dello spazio di Hilbert"
                        "\n\t2.\"nome file rho iniziale\"\n\t3.\"nome file "
                        "H\"\n\t4. numero di operatori d'errore\n\t5.\""
                        "nome file operatori d'errore\"\n\t6. \"nome file "
                        "ratei errore\"\n\t7. numero osservabili\n\t8. "
                        "\"nome file osservabili\"\n\t9. tempo finale della "
                        "simulazione\n\t10. numero di intervalli temporali "
                        "mesh\n\t11. frequenza con cui stampare la rho del "
                        "sistema\n\t12. nome file in cui stampare "
                        "l'avanzamento della simulazione\n\n");
        return 1;
    }

/*
 *  inizializziamo i valori passati in input
 */
    uint32_t dim   = atoi(argv[0]);
    uint32_t n_err = atoi(argv[3]);
    uint32_t n_oss = atoi(argv[6]);
    
    double t_fin   = atof(argv[8]);
    double N       = atof(argv[9]);
    double N2      = atof(argv[10]);

    char *file_rho_in = argv[1];
    char *file_H      = argv[2];
    char *file_err    = argv[4];
    char *file_ratei  = argv[5];
    char *file_oss    = argv[7];
    char *file_av     = argv[11];

/*
 *  Allochiamo le strutture dati necessarie per l'evoluzione temporale del
 *  sistema
 */
    parametri_simulazione *param = 
        alloc_parametri_simulazione(t_fin, (uint32_t)N, (uint32_t)N2,
                                    n_oss, dim, file_av);
    spazio_lavoro_forzante *slf =
        alloc_spazio_lavoro_forzante(dim, n_err);

/*
 *  Verifichiamo l'esito delle allocazioni
 */
    if (param == NULL || slf == NULL)
    {
        free_parametri_simulazione(&param);
        free_spazio_lavoro_forzante(&slf);
        return 2;
    }
/*
 *  Leggiamo i file passati in input
 */
    uint8_t ris = 
        leggi_input(param, slf, file_rho_in, file_H, file_err, file_ratei, file_oss);
/*
 *  Verifichiamo l'esito della lettura dei file
 */
    if (ris > 0)
    {
        free_parametri_simulazione(&param);
        free_spazio_lavoro_forzante(&slf);
        return 3;
    }

/*
 *  Eseguiamo l'integrazione numerica del sistema di ode del problema
 */
    simula(param, slf);

/*
 *  Liberiamo la memoria precedentemente allocata
 */
    free_parametri_simulazione(&param);
    free_spazio_lavoro_forzante(&slf);

    return 0;
}