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
                                            n: int,                     \
                                            layout: list[int],          \
                                            backend,                    \
                                            fold_factors: list[int],    \
                                            out):                       \
#   Funzione che calcola i valori di aspettazione delle varie potenze
#   dell'Hamiltoniana necessarie per la costruzione di una Krylov-like 
#   subspace expansion. Ogni valore di aspettazione è calcolato per diverse
#   intensità d'errore, in modo da poter eseguire una ZNE su tali valori
    
    #   Calcoliamo la potenza massima di cui calcolare il valore di aspettazione
    pow_max = (2 * n) + 1

#   Calcoliamo il primo grouping degli elementi che commutano di Ham
    # oss = Ham.tot

#   Output di avanzamento
    print("Inizio stima valore di aspettazione della 1-potenza dell'Hamiltoniana")

#   Stimiamo il valore di aspettazione di H
    # stima_aspettazione_ZNE(circuit, oss, estimator, fold_factors, layout, backend, out)

#   Cicliamo su tutte le potenze di cui calcolare i valori d'aspettazione
    for i in range(2, pow_max + 1):

#       Calcoliamo la potenza attuale dell'Hamiltoniana suddivisa in gruppi di
#       operatori che commutano tra loro
        # oss = (oss @ Ham.tot).simplify()

#       Output di avanzamento
        print("Inizio stima valore di aspettazione della %d-potenza dell'Hamiltoniana" % i)

#       calcoliamo il valore di aspettazione dell'osservabile attuale
        # stima_aspettazione_ZNE(circuit, oss, estimator, fold_factors, layout, backend, out)
    
#   end i

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
    va = np.zeros(pow_max + 1)

#   Valore di aspettazione dell'identità
    va[0] = 1

#   Calcoliamo il primo grouping degli elementi che commutano di Ham
    oss = Ham.tot

#   Scaling factors per la ZNE
    fold_factors = [1, 1.5, 2, 2.5, 3, 4, 5]

#   DA RIPENSARE:: COME INCLUDERE LA ZNE
#   Stimiamo il valore di aspettazione di H
    va[1] = stima_aspettazione_ZNE(circuit, oss, estimator, fold_factors, layout, backend, out)

#   Cicliamo su tutte le potenze di cui calcolare i valori d'aspettazione
    for i in range(2, pow_max + 1):

#       Calcoliamo la potenza attuale dell'Hamiltoniana suddivisa in gruppi di
#       operatori che commutano tra loro
        oss = (oss @ Ham.tot).simplify()

#       calcoliamo il valore di aspettazione dell'osservabile attuale
        va[i] += stima_aspettazione_ZNE(circuit, oss, estimator, fold_factors, layout, backend, out)
    
#   end i

#   Aggiorniamo l'output su file ed eseguiamo un flush di questo
    out.write("\n----------------------------------------\n")
    out.write("Terminata costruzione SbE")
    out.write("\n----------------------------------------\n")
    out.write("\n\n")
    out.flush()

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