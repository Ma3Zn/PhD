from povm_toolbox.library import *
from qiskit.quantum_info import *
from qiskit.primitives import *
from qiskit.providers import *
from qiskit import *
from HamOp import *
from Circ import *
from math import *
from SbE import *

# ----------------------- Fixed System parameters --------------------------- #

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

# --------------------------------------------------------------------------- #

# Function that generate a SparsePauliOp representation of the Hamiltonian of 
# the system
def build_hamiltonian(n_spin: int) -> HamOp:
    return crea_hamiltonianaOp(n_spin, chiuso)


# Function that compute a numerical value for the Hamiltoniana of the system
def evaluate_hamiltonian(H: HamOp):
    return H.valuta_tot(J, drive, lam, 1/lam)

# Function that generate the circuit to perform the proper digitalization of
# the annealing process
def trotterized_annealing(n_spin: int, n_tstep: int) -> QuantumCircuit:
    
    init = crea_circuito_inizializzazione(n_spin)

    digitalized_annealing = crea_circuito_annealing(J, n_spin, n_tstep, \
                                                    drive, lam, drive_sp, lam)
    
    return init.compose(digitalized_annealing)


# Function that compute an estimate for the ground state energy by a Subspace 
# Expansion
def energy_estimate_SbE(circ: QuantumCircuit,       \
                        H: HamOp,                   \
                        POVM_result: POVMPubResult, \
                        max_pow: int)               \
                        -> float:
    return stima_energia_SbE(circ, H, POVM_result, max_pow)