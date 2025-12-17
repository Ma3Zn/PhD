from qiskit_ibm_runtime import EstimatorV2 as Estimator
from qiskit_ibm_runtime import SamplerV2

import scipy.linalg as ln

import mapomatic as mm

from Ring import *

from qiskit_ibm_runtime.fake_provider import FakeTorino, FakeBrisbane

from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

# Ring sites
n_spin = 12

# Creation of the system Hamiltonian as SparsePauliOp
H = build_hamiltonian(n_spin)

# Number of trotter step -- each trotter step is 15-CNOT deep
n_tstep = 4

# Creation of the annealing circuit
qc = trotterized_annealing(n_spin, n_tstep)

# Selection of the fake backend
fake_backend = FakeTorino()

# # TEST QPU
# from qiskit_ibm_runtime import QiskitRuntimeService

# # Mio Token
# service = QiskitRuntimeService(
#     channel='ibm_quantum',
#     instance='ibm-q/open/main',
#     token='feb77ce69166caff54852530c1634b44e307ee4fd5af82cbea78c06774b07ed575a010878366428456e3333366c38aed28bc7a6dd61279eca7febcebc353f3b3'
# )

# # Selezione manuale backend
# fake_backend = service.backend("ibm_brisbane")

# Transpiling the circuit for the fake backend
trans_qc = transpile(qc, fake_backend, optimization_level=3)

# Reducing the size of the circuit
small_qc = mm.deflate_circuit(trans_qc)

# Finding the matchin subgraph for the fake backend
layouts = mm.matching_layouts(small_qc, fake_backend)

# Evaluating the cost of each layout
scores = mm.evaluate_layouts(small_qc, layouts, fake_backend)

print(scores)

# Creation of the standard pass menager for the transpiling pipeline with the 
# best possible initial layout
pm = generate_preset_pass_manager(optimization_level=3, \
                                  backend=fake_backend, \
                                  initial_layout=scores[0][0])

# Transpile the circuit to an "Instruction Set Architecture" (ISA) circuit.
# Note: the transpiler automatically adds "ancilla" qubits to make the
# transpiled circuit match the size of the FakeSherbrooke backend.
isa_circ = pm.run(small_qc)

# Setting the maximum power of the Hamiltonian to use
max_pow = 1

# Allocating the estiamtor to use for the Subpspace Expansion
estimator = Estimator(mode=fake_backend)

# # EXACT Estimator
# estimator = Estimator()

# Retriving an energy estimation throught digitalized annealing + 
# subspace expansion
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