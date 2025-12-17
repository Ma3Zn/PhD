from qiskit import *
from qiskit.quantum_info import *
from qiskit_ibm_runtime import Estimator

from HamOp import *

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


def genera_matrici_SbE(circuit: QuantumCircuit,         \
                       Ham: HamOp,                      \
                       estimator: Estimator,            \
                       n: int):
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

#   Applichiamo il layout opportuno all'osservabile [NON serve]
    isa_oss = oss.apply_layout(circuit.layout)

#   Stimiamo il valore di aspettazione di H
    va[1] = stima_aspettazione(circuit, isa_oss, estimator)

#   Cicliamo su tutte le potenze di cui calcolare i valori d'aspettazione
    for i in range(2, pow_max + 1):

#       Calcoliamo la potenza attuale dell'Hamiltoniana suddivisa in gruppi di
#       operatori che commutano tra loro
        oss = (oss @ Ham.tot).simplify()

#       Applichiamo il layout opportuno all'osservabile [NON serve]
        isa_oss = oss.apply_layout(circuit.layout)

#       calcoliamo il valore di aspettazione dell'osservabile attuale
        va[i] += stima_aspettazione(circuit, isa_oss, estimator)
    
#       end k
#   end i

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
                        n: int)                 \
                        -> float:
#   Funzione che data un'Hamiltoniana ed il numero di potenze con cui generare
#   lo spazio di Krylov, esegue una Subspace Expansion per recuperare una stima
#   migliorata del valore energetico del ground state ottenuto tramite
#   il circuito passato

#   Generiamo le matrici del problmea agli autovalori generalizzato che nasce
#   dalla Subspace Expansion (pesudocodice)
    H, S = genera_matrici_SbE(circuit, Ham, estimator, n)

#   Risolviamo il problema agli autovalori generalizzato e ritorniamo
#   l'autovalore più piccolo in modulo
    return ln.eigh(H, S, eigvals_only=True, subset_by_index=[0, 0])