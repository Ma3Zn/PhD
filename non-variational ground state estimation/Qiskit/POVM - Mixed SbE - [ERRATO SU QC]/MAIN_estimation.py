from qiskit.primitives import Estimator

import scipy.linalg as ln

from povm_toolbox.library import *

from Ring import *

# Ring sites
n_spin = 12

# Creation of the system Hamiltonian as SparsePauliOp
H = build_hamiltonian(n_spin)

# Number of trotter step -- each trotter step is 15-CNOT deep
n_tstep = 7

# Creation of the annealing circuit
circ = trotterized_annealing(n_spin, n_tstep)

# Setting the maximum power of the Hamiltonian to use
max_pow = 1

# Excitaion operators
(op, n_op) = build_exitation_operators(n_spin)

# Allocating the POVMscheme to use for the expectation value estimation
POVM_scheme = ClassicalShadows(num_qubits=n_spin, measurement_twirl=False, shot_repetitions=1)
POVM_scheme.measurement_circuit.draw("mpl")

# Retriving an energy estimation throught digitalized annealing + subspace expansion
estimate = energy_estimate_SbE(circ, H, op, POVM_scheme, max_pow, n_op)

# Computing a numerical value of the system Hamiltonian
H_num = evaluate_hamiltonian(H)

# Numerical diagonalization of the Hamiltonian
e, V = ln.eigh(np.real(H_num), subset_by_index=[0, 0])

# Recuperiamo l'indice del groung state
id_gs = np.argmin(e)

# Recuperiamo il ground state e la relativa energia
gs     = V[:,id_gs]
energy = e[id_gs]

err = abs(energy - estimate) / abs(energy)
print([energy, estimate, err])