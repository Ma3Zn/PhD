from qiskit.primitives import BackendEstimatorV2 as Estimator
from qiskit.quantum_info import *
from qiskit import *
from CustomBackend import *
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
# #   Drive lineare
#    ris = param * t

#   Drive sinusoidale
    ris = sin((param * t) * pi/2) ** (0.9)

    return ris

# Funzione di drive per l'annealing del termine a singolo corpo
def drive_sp(param, n, t):
# #   Potenza
#    ris = 0

# #    if ((param * t)**n < 1):
#        ris = 1 - (param * t)**n

#   Tangente
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


# Function that compute the expectation values needed for the construction of 
# a krylov-like subspace expansion at different noise-strenght levels 
# (for a successive ZNE). It writes such expectation values on file.
def compute_expextation_values_SbE_for_ZNE(circ: QuantumCircuit,        \
                                            H: HamOp,                   \
                                            estimator: Estimator,       \
                                            min_pow: int,               \
                                            max_pow: int,               \
                                            layout: list[int],          \
                                            backend,                    \
                                            folding_factors: list[int], \
                                            out_mean, out_var, out_ff): \
    calcola_valori_aspettazione_SbE_per_ZNE(circ, H, estimator, min_pow, max_pow, layout, backend, folding_factors, out_mean, out_var, out_ff)

# Function that compute an estimate for the ground state energy by a Subspace 
# Expansion
def energy_estimate_SbE(circ: QuantumCircuit,   \
                        H: HamOp,               \
                        estimator: Estimator,   \
                        max_pow: int,           \
                        layout: list[int],      \
                        backend,                \
                        out)                    \
                        -> float:
    return stima_energia_SbE(circ, H, estimator, max_pow, layout, backend, out)
