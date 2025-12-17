from qiskit import *
from qiskit.quantum_info import *
from qiskit_ibm_runtime import EstimatorV2 as Estimator

from HamOp import *

from ZNE import stima_aspettazione_ZNE

import scipy.linalg as ln
import numpy as np

def stima_aspettazione(circuit: QuantumCircuit,     \
                       oss: SparsePauliOp,          \
                       estimator: Estimator)        \
                       -> float:
#   Funzione che dato un circuito, un'osservabile ed un estimator, stima il 
#   valore di aspettazione  di questa

#   Creazione ed esecuzione del job
    job = estimator.run([(circuit, oss)])
    pub_result = job.result()[0]

#   Ritorniamo il valore d'aspettazione stimato
    return pub_result.data.evs


def calcola_valori_aspettazione_SbE_per_ZNE(circuit: QuantumCircuit,    \
                                            Ham: HamOp,                 \
                                            estimator: Estimator,       \
                                            pow_min: int,               \
                                            pow_max: int,               \
                                            layout: list[int],          \
                                            backend,                    \
                                            fold_factors: list[int],    \
                                            out_mean, out_var, out_ff): \
#   Funzione che calcola i valori di aspettazione delle varie potenze
#   dell'Hamiltoniana necessarie per la costruzione di una Krylov-like 
#   subspace expansion. Ogni valore di aspettazione è calcolato per diverse
#   intensità d'errore, in modo da poter eseguire una ZNE su tali valori

#   Controlliamo se la potenza minima di cui calcolare il valore 
#   d'aspettaione è l'identità o meno
    if (pow_min <= 0):

#       Stampiamo un messaggio di warning all'utente
        print("\nWARNING:: potenza minima richiesta troppo piccola.\n "
              "Il valore di aspettazione dell'identita' e' 1.\n")

#       Settiamo il valore minimo della potenza di cui calcolare il valore 
#       d'aspettazione ad 1
        pow_min = 1

#   lista osservabili da calcolare
    lst_oss = []

#   Calcoliamo la prima potenza dell'osservabile
    oss = Ham.tot

#   Controlliamo se dobbiamo includerla o meno
    if (pow_min == 1):

#       Aggiungiamo l'osservabile alla lista delle osservabili da calcolare
        lst_oss.append(oss)
        
#       Output di avanzamento
        print("Stima valore di aspettazione della 1-potenza dell'Hamiltoniana")

#   Cicliamo su tutte le potenze di cui calcolare i valori d'aspettazione
    for i in range(2, pow_max + 1):

#       Calcoliamo la potenza attuale dell'Hamiltoniana suddivisa in gruppi di
#       operatori che commutano tra loro
        oss = (oss @ Ham.tot).simplify()

#       Controlliamo se dobbiamo stimare o meno il valore di aspettazione della potenza 
#       attuale
        if (i >= pow_min):

#           Aggiungiamo l'osservabile alla lista delle osservabili da calcolare
            lst_oss.append(oss)

#           Output di avanzamento
            print("Stima valore di aspettazione della %d-potenza dell'Hamiltoniana" % i)
    
#   end i

#   Siccome H commuta sempre con tutte le sue potenze (questo è vero per un 
#   qualsiasi operatore) calcoliamo i vari valori di aspettazione in un'unica
#   chiamata all'estimator per tutti i diversi folding_factors
    stima_aspettazione_ZNE(circuit, lst_oss, estimator, fold_factors, layout, backend, out_mean, out_var, out_ff)

def genera_matrici_SbE(circuit: QuantumCircuit, \
                       Ham: HamOp,              \
                       estimator: Estimator,    \
                       n: int,                  \
                       layout: list[int],       \
                       backend,                 \
                       out):
#   Funzione che passato stimatore, circuito, Hamiltoninaa e numero di potenze
#   di cui calcolare i valori di aspettazione genera le matrici relative al
#   problema agli autovalori generalizzato generato dalla Subspace Expansion

#   Calcoliamo le dimensioni delle matrici
    dim = n+1

#   Allochiamo le matrici
    H = np.zeros((dim, dim))
    S = np.zeros((dim, dim))

#   Calcoliamo la potenza massima di cui calcolare il valore di aspettazione
    pow_max = (2 * n) + 1

#   Allochiamo un vettore che conterrà tutti i valori di aspettazione calcolati
    va = []

#   Valore di aspettazione dell'identità
    va.append(1)

#   lista osservabili da calcolare
    lst_oss = []

#   Calcoliamo il primo grouping degli elementi che commutano di Ham
    oss = Ham.tot

#   Aggiungiamo l'osservabile attuale alla lista delle osservabili da calcolare
    lst_oss.append(oss)

#   Cicliamo su tutte le potenze di cui calcolare i valori d'aspettazione
    for i in range(2, pow_max + 1):

#       Calcoliamo la potenza attuale dell'Hamiltoniana suddivisa in gruppi di
#       operatori che commutano tra loro
        oss = (oss @ Ham.tot).simplify()

#       Aggiungiamo l'osservabile attuale alla lista delle osservabili da calcolare
        lst_oss.append(oss)
    
#   end i

#   Siccome H commuta sempre con tutte le sue potenze (questo è vero per un 
#   qualsiasi operatore) calcoliamo i vari valori di aspettazione in un'unica
#   chiamata all'estimator
    tmp = stima_aspettazione(circuit, lst_oss, estimator)
    va = [*va, *tmp]

#   Con i valori di aspettazione ottenuti costruimao le matrici del problema
#   agli autovalori generalizzato
    for i in range(0,dim):
        for j in range(0,dim):
#           Recuperiamo l'indice della potenza attuale [per H --> @-1 per S]
            pow_idx = 1 + i + j

#           Inseriamo l'opportuno valore di aspettazione in H
            H[i,j] = va[pow_idx]

#           Inseriamo l'opportuno valore di aspettazione in S
            S[i,j] = va[i + j]
#       end j
#   end i

    return H, S

def stima_energia_SbE(circuit: QuantumCircuit,  \
                        Ham: HamOp,             \
                        estimator: Estimator,   \
                        n: int,                 \
                        layout: list[int],      \
                        backend,                \
                        out)                    \
                        -> float:
#   Funzione che data un'Hamiltoniana ed il numero di potenze con cui generare
#   lo spazio di Krylov, esegue una Subspace Expansion per recuperare una stima
#   migliorata del valore energetico del ground state ottenuto tramite
#   il circuito passato

#   Generiamo le matrici del problmea agli autovalori generalizzato che nasce
#   dalla Subspace Expansion (pesudocodice)
    H, S = genera_matrici_SbE(circuit, Ham, estimator, n, layout, backend, out)

#   Risolviamo il problema agli autovalori generalizzato e ritorniamo
#   l'autovalore più piccolo in modulo
    return ln.eigh(H, S, eigvals_only=True, subset_by_index=[0, 0])