from qiskit_ibm_runtime import EstimatorV2 as Estimator
from qiskit_ibm_runtime import SamplerV2

import scipy.linalg as ln

from Ring import *

from qiskit_ibm_runtime.fake_provider import FakeBrisbane

fake_backend = FakeBrisbane()

# Fixing the qubits to be used on the QPU 
layout = [37, 38, 39, 40 ,41, 53, 60, 59, 58, 57, 56, 52]

# Fanne uno per ognuno degli hardware in modo da utilizzare
# sempre il migliore possibile [TODO]

# Puoi anche investigare il plug_in "mapomatic" che vuole
# automatizzare questo step. [TODO]
# https://github.com/qiskit-community/mapomatic

# Ring sites
n_spin = 12

# Creation of the system Hamiltonian as SparsePauliOp
H = build_hamiltonian(n_spin)

# Number of trotter step -- each trotter step is 15-CNOT deep
n_tstep = 4

# Creation of the annealing circuit
qc = trotterized_annealing(n_spin, n_tstep)

from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
pm = generate_preset_pass_manager(optimization_level=2, backend=fake_backend, initial_layout=layout)

# Transpile the circuit to an "Instruction Set Architecture" (ISA) circuit.
# Note: the transpiler automatically adds "ancilla" qubits to make the transpiled
# circuit match the size of the FakeSherbrooke backend.
isa_circ = pm.run(qc)

# circ.draw()

# Setting the maximum power of the Hamiltonian to use
max_pow = 1

# # Allocating the estiamtor to use for the Subpspace Expansion -- EXACT Estimator
# estimator = Estimator()
estimator = Estimator(mode=fake_backend)

# Retriving an energy estimation throught digitalized annealing + subspace expansion
estimate = energy_estimate_SbE(isa_circ, H, estimator, max_pow)

# Computing a numerical value of the system Hamiltonian [TBD]
H_num = evaluate_hamiltonian(H)

# Numerical diagonalization of the Hamiltonian [TBD]
e, V = ln.eigh(np.real(H_num), subset_by_index=[0, 0])

# Recuperiamo l'indice del groung state
id_gs = np.argmin(e)

# Recuperiamo il ground state e la relativa energia
gs     = V[:,id_gs]
energy = e[id_gs]

err = abs(energy - estimate) / abs(energy)
print([energy, estimate, err])