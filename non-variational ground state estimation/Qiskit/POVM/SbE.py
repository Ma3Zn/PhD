from qiskit import *
from qiskit.primitives import *
from qiskit.quantum_info import *
from povm_toolbox.library import *
from povm_toolbox.sampler import *
from povm_toolbox.post_processor import POVMPostProcessor

from HamOp import *
from POVM import *

import scipy.linalg as ln
import numpy as np

def genera_matrici_SbE(circuit: QuantumCircuit,         \
                       Ham: HamOp,                      \
                       POVM_result: POVMPubResult,      \
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

#   Allochiamo un vettore che conterrà tutte le deviazioni standard delle stime
#   dei valori di aspettazioni prodotte dalle POVM
    std = np.zeros(pow_max + 1)

#   DEBUG:: variabike che contine i valori di aspettazione esatti
    va_esatti = np.zeros(pow_max + 1)

#   Allochiamo il post processor per la stima dei valori di aspettazione
#   tramite POVM
    post_processor = alloca_POVM_post_processor(POVM_result)

#   Valore di aspettazione dell'identità
    va[0] = 1
    std[0] = 0
    va_esatti[0] = 1 # DEBUG

#   Calcoliamo il primo grouping degli elementi che commutano di Ham
    oss = Ham.tot
    
#   Stimiamo il valore di aspettazione di H. Al momento non ci interessiamo
#   delle std_dev stimate dal metodo per i vari valori di aspettazione
    va[1], std[1] = post_processor.get_expectation_value(oss)
    va_esatti[1] = Statevector(circuit).expectation_value(oss)

#   Cicliamo su tutte le potenze di cui calcolare i valori d'aspettazione
    for i in range(2, pow_max + 1):

#       Calcoliamo la potenza attuale dell'Hamiltoniana suddivisa in gruppi di
#       operatori che commutano tra loro
        oss = (oss @ Ham.tot).simplify()

#       calcoliamo il valore di aspettazione dell'osservabile attuale
        va[i], std[i] = post_processor.get_expectation_value(oss)
        va_esatti[i] = Statevector(circuit).expectation_value(oss) # DEBUG
#   end i

#   DEBUG:: Stampiamo l'errore relativo sulla stima dell'osservabile
    for i in range(0, pow_max + 1):
        print(f"errore: {abs(va[i] - va_esatti[i]):>14.6f}")
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

def stima_energia_SbE(circuit: QuantumCircuit,      \
                        Ham: HamOp,                 \
                        POVM: POVMImplementation,   \
                        n: int)                     \
                        -> float:
#   Funzione che data un'Hamiltoniana ed il numero di potenze con cui generare
#   lo spazio di Krylov, esegue una Subspace Expansion per recuperare una stima
#   migliorata del valore energetico del ground state ottenuto tramite
#   il circuito passato

#   Eseguiamo le POVM
    POVM_result = esegui_POVM(circuit, POVM)

#   Generiamo le matrici del problmea agli autovalori generalizzato che nasce
#   dalla Subspace Expansion (pesudocodice)
    H, S = genera_matrici_SbE(circuit, Ham, POVM_result, n)

#   Risolviamo il problema agli autovalori generalizzato e ritorniamo
#   l'autovalore più piccolo in modulo
    return ln.eigh(H, S, eigvals_only=True, subset_by_index=[0, 0])