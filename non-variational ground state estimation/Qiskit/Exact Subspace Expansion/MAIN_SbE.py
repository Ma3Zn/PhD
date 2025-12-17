from HamOp import *
from Circ import *
from SbE import *

import numpy as np
from numpy import linalg

import scipy.linalg as ln

import matplotlib
import matplotlib.pyplot as plt

from math import *

from qiskit import *
from qiskit_aer import QasmSimulator, UnitarySimulator, StatevectorSimulator
from qiskit.quantum_info import Pauli

from qiskit.primitives import Estimator

from qiskit.visualization import plot_bloch_multivector

# Usiamo il BaseEstimatorV1 di qiskit
estimator = Estimator()

# Numero spin anello
n_spin = 4

# Periodicià anello
chiuso = True

# Costanti d'interazione [in unità di J]
J = [1,1,1]

# Parametro per la regolazione della velocità del processo di annealing
lam = 0.1

# Funzione di drive per l'annealing
def drive(param, t):
#   Drive lineare
    ris = sin((param * t) * pi/2) ** (0.9)

    return ris

# Funzione di drive per l'annealing del termine a singolo corpo
def drive_sp(param, n, t):
    ris = tan(1) - tan((param * t) ** n)

    return ris

# Creiamo l'Hamiltoniana del sistema
H = crea_hamiltonianaOp(n_spin, chiuso)

# Creiamo il circuito per la preparazione dello stato iniziale
init = crea_circuito_inizializzazione(n_spin)
#init.draw()

# Numero di trotter step in cui digitalizzare il processo di annealing
n_tstep = 4

# Creaimo il circuito per l'esecuzione dell'annealing digitalizzato
digitalized_annealing = crea_circuito_annealing(J, n_spin, n_tstep, drive, lam, drive_sp, lam)
#digitalized_annealing.draw()

# Uniamo i circuiti derivanti dall'inizializzazione e dal processo di annealing
final_circuit = init.compose(digitalized_annealing)
# final_circuit.draw()

# Allochiamo il backend opportuno, i.e. lo statevector_simulator
backend = StatevectorSimulator()

# Eseguiamo la simulazione e recuperiamo lo statevector del nostro sistema
# approssimazione del ground state
psi_fin = backend.run(final_circuit).result().get_statevector()

# Calcolo energia psi finale dopo annealing
e_annealing = stima_energia_SbE(final_circuit, H, estimator, 0)

# ATTENZIONE:: per Qiskit |0> == |up> e |1> == |down>
# plot_bloch_multivector(psi_fin)

# Valutazione Hamiltoniana di Heisenberg del sistema
H_fin = H.valuta_tot(J, drive, lam, 1/lam)

# Diagonalizzazione Hamiltoninana di Heisenberg
e, V = ln.eigh(np.real(H_fin), subset_by_index=[0, 0])

# Recuperiamo l'indice del groung state
id_gs = np.argmin(e)

# Recuperiamo il ground state e la relativa energia
gs   = V[:,id_gs]
e_gs = e[id_gs]

# Numero potenze dell'Hamiltoniana da utilizzare per la costruzione del SbE
n_potenze = 2

e_SbE = stima_energia_SbE(final_circuit, H, estimator, n_potenze)
print(e_SbE)

# Calcolo errore relativo per l'energia dopo la sola fase di annealing
err_e_annealing = abs(e_gs - e_annealing) / abs(e_gs)

# Calcolo errore relativo per l'energia dopo la SbE
err_e_SbE = abs(e_gs - e_SbE) / abs(e_gs)

# Miglioramento relativo dell'errore con l'ausilio della SbE
gain_SbE = abs(err_e_annealing - err_e_SbE) / abs(err_e_annealing)

# Stampiamo i risultati
print([e_gs, e_SbE[0], e_annealing[0]])
print([err_e_annealing[0], err_e_SbE[0], gain_SbE[0]])