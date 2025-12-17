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
import math

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

def cpx_clean(num: complex, tol: float = 1e-14) -> complex:
#   Funzione che rimuove eventuale sporcizia numerica da un numero complesso

#   Controlliamo la parte reale
    if (abs(num.real) < tol):
        ris = num * 0 + num.imag * 1j

#   Controlliamo la parte immaginaria
    if (abs(num.imag) < tol):
        ris = num.real + num * 0j

#   Ritorniamo il numero "pulito"
    return ris

#-----------------------------------------------------------------------------#

def crea_operatori_eccitazione_SbE(n_spin: int, loc: bool = False) \
    -> tuple[list[SparsePauliOp], int]:
#   Funzione che dato il numero di sisti nell'anello genera gli operatori di
#   eccitazione da utilizzare per la costruzione della subspace expansion

#   Contatore operatori
    n = 1

#   inizializzazione lista operatori
    op = []

#   Inseriamo il primo operatorore --> l'identità
#   Risulta necessario includere anche questo operatore per mantenere anche il
#   sottospazio generato dalle potenze di H desiderate
    op.append(SparsePauliOp(''.join(['I']*n_spin)))

#   Operatori d'eccitazione::   Sz locale
    if loc:
#       generazione stringa operatore di pauli
        cod = ['I'] * n_spin

#       Inserimento in posizione opportuna di Sz
        cod[0] = 'Z'
        
#       Generazione dell'opportuno operatore di Pauli
        op.append(SparsePauliOp(''.join(cod), 1/2))

#       aggiornamento numero di operatori inseriti
        n = n + 1

#   Operatori d'eccitazione::   combinazione Sz [TODO]
    else:
#       Calcoliamo tutti i valori di k ammissibili
        k_vals = np.linspace(0, 2 * math.pi, n_spin + 1)
        k_vals = k_vals[1:]

#       Cicliamo sui valori di k
        for k in k_vals:
#           Aggiorniamo il contatore degli operatore generati
            n = n + 1

#           Inizializziamo la lista delle stringhe degli operatori
            lista_stringhe_op = []

#           Inizializziamo il vettore dei pesi
            pesi = np.zeros(n_spin, dtype=complex)

#           Cicliamo sui vari siti
            for i in range(0,n_spin):

#               generiamo la stringa iniziale
                cod = ['I'] * n_spin

#               Inseriamo l'operatore nel sito opportuno
                cod[i] = 'Z'

#               Appendiamo l'operatore attuale alla lista degli operatori
                lista_stringhe_op.append(''.join(cod))

#               Inseriamo il peso dell'operatore nel vettore dei pesi
                pesi[i] = cpx_clean(np.exp(1j * k * i)/2)
#           end i

#           Generiamo l'opportuno operatore
            op.append(SparsePauliOp(lista_stringhe_op, pesi))
#       end k

#   Ritorniamo gli operatori generati e il loro numero
    return (op, n)

#-----------------------------------------------------------------------------#

def costruisci_operatori_SbE(Ham: HamOp,                \
                             op: list[SparsePauliOp],   \
                             pow_max: int,                    \
                             n_op: int)                 \
                             -> tuple[list[SparsePauliOp], int]:
#   Funzione che costruisce il set totale di operatori per la SbE combinando
#   potenze dell'Hamiltonina e gli operatori passati

#   Allocazione lista operatori
    operator_set = []

#   Costruiamo le potenze dell'Hamiltoniana [TODO]
    Ham_pow = costruisci_potenze_Ham(Ham, pow_max)

#   Cicliamo sulle potenze dell'Hamiltoniana
    for i in range(0, pow_max):
#       Cicliamo sugli operatori
        for j in range(0, n_op):

#           Aggiungiamo l'operatore desiderato alla lista
            operator_set.append(Ham_pow[i] @ op[j])

        # end j
    #end i

    return (operator_set, len(operator_set))

#-----------------------------------------------------------------------------#

def genera_matrici_SbE(circuit: QuantumCircuit,         \
                       Ham: HamOp,                      \
                       op: list[SparsePauliOp],         \
                       POVM_result: POVMPubResult,      \
                       n: int,                          \
                       n_op: int):
#   Funzione che passato stimatore, circuito, Hamiltoninaa, numero di potenze
#   ed operatori di eccitazione di cui calcolare i valori di aspettazione
#   genera le matrici relative al problema agli autovalori generalizzato
#   generato dalla Subspace Expansion

#   Calcoliamo le dimensioni delle matrici
    dim = (n + 1) * n_op

#   Allochiamo le matrici
    H = np.zeros((dim, dim))
    S = np.zeros((dim, dim))

#   Costruiamo gli operatori opportuni alla relizzazione della SbE
    (operator_set, n_op_tot) = costruisci_operatori_SbE(Ham, op, n, n_op)

#   DEBUG:: variabile che conterrà le deviazioni standard stimate per i valori
#           di aspettazioni approssimati tramite POVM
    std = [0, 0]

#   DEBUG:: variabile che contine i valori di aspettazione esatti
    va_esatti = [0, 0]

#   Allochiamo il post processor per la stima dei valori di aspettazione
#   tramite POVM
    post_processor = alloca_POVM_post_processor(POVM_result)

#   DEBUG:: Output su precisione dei valori di aspettaione
    print("(H)   Estimated value   Estimated std   Actual error       \
           (s)   Estimated value   Estimated std   Actual error")

#   Cicliamo su tutti gli operatori della SbE per custruire le matrici del
#   problema generalizzato agli autovalori
    for i in range(0, n_op_tot):
        for j in range(0, n_op_tot):

#           Costruzione di H
            op_H = operator_set[i].adjoint() @ Ham.tot @ operator_set[j]

#           DEBUG:: rapprsentazione matriciale operatore
            tmp_H = abs(op_H.to_matrix() - op_H.adjoint().to_matrix())

#           calcoliamo il valore di aspettazione dell'osservabile attuale
            H[i,j], std[0] = post_processor.get_expectation_value(op_H)

#           DEBUG:: Calcolo del valore di aspettazione esatto
            va_esatti[0] = Statevector(circuit).expectation_value(op_H)

#           Costruzione di S
            op_S = operator_set[i].adjoint() @ operator_set[j]

#           DEBUG:: rapprsentazione matriciale operatore
            tmp_s = op_S.to_matrix()

#           calcoliamo il valore di aspettazione dell'osservabile attuale
            S[i,j], std[1] = post_processor.get_expectation_value(op_S)

#           DEBUG:: Calcolo del valore di aspettazione esatto
            va_esatti[1] = Statevector(circuit).expectation_value(op_S)

# #           DEUBG:: Outpunt su precisione dei valori di aspettazione
#             print(f"{H[i,j]:>17.6f} {std[0]:>15.6f} \
#                     {abs(H[i,j] - va_esatti[0]):>14.6f}\t\
#                     {S[i,j]:>17.6f} {std[1]:>15.6f} \
#                     {abs(S[i,j] - va_esatti[1]):>14.6f}")
            
#           DEBUG:: Utilizziamo i valori di aspettazione esatti per la SbE
            H[i,j] = va_esatti[0]
            S[i,j] = va_esatti[1]

#       end j
#   end i

    return H, S

#-----------------------------------------------------------------------------#

def stima_energia_SbE(circuit: QuantumCircuit,      \
                        Ham: HamOp,                 \
                        op: list[SparsePauliOp],    \
                        POVM: POVMImplementation,   \
                        max_pow: int,               \
                        n_op: int)                  \
                        -> float:
#   Funzione che data un'Hamiltoniana ed il numero di potenze con cui generare
#   lo spazio di Krylov, esegue una Subspace Expansion per recuperare una stima
#   migliorata del valore energetico del ground state ottenuto tramite
#   il circuito passato

#   Eseguiamo le POVM
    POVM_result = esegui_POVM(circuit, POVM)

#   Generiamo le matrici del problmea agli autovalori generalizzato che nasce
#   dalla Subspace Expansion (pesudocodice)
    H, S = genera_matrici_SbE(circuit, Ham, op, POVM_result, max_pow, n_op)

#   Risolviamo il problema agli autovalori generalizzato e ritorniamo
#   l'autovalore più piccolo in modulo
    return ln.eigh(H, S, eigvals_only=True, subset_by_index=[0, 0])

#-----------------------------------------------------------------------------#