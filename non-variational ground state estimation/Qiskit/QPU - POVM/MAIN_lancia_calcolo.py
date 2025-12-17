from qiskit.primitives import Estimator

import scipy.linalg as ln

from povm_toolbox.library import *

from Ring import *

from POVM import lancia_POVM

from qiskit_ibm_runtime import QiskitRuntimeService

from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

# Script che lancia su QPU IBM l'esecuzione del circuito con POVM

# # Mio Token
# service = QiskitRuntimeService(
#     channel='ibm_quantum',
#     instance='ibm-q/open/main',
#     token='feb77ce69166caff54852530c1634b44e307ee4fd5af82cbea78c06774b07ed575a010878366428456e3333366c38aed28bc7a6dd61279eca7febcebc353f3b3'
# )

# Token Alessandro
service = QiskitRuntimeService(
    channel='ibm_quantum',
    instance='ibm-q/open/main',
    token='a67cbd727256bbab2306b1857d9a0fcbc82dc599a32a20157a0be3e76894a4a8860a4c69252121fbb7f03d9caca6f126661885cddbc21c0bc39a5eb3cf62220d'
)

# # Baackend meno utilizzato
# backend = service.least_busy(operational=True, min_num_qubits=12)
# print(
#     f"Name: {backend.name}\n"
#     f"Version: {backend.version}\n"
#     f"No. of qubits: {backend.num_qubits}\n"
# )

# Selezione manuale backend
backend = service.backend("ibm_brisbane")

# Ring sites
n_spin = 12

# Creation of the system Hamiltonian as SparsePauliOp
H = build_hamiltonian(n_spin)

# Number of trotter step -- each trotter step is 15-CNOT deep
n_tstep = 4

# Creation of the annealing circuit
circ = trotterized_annealing(n_spin, n_tstep)

# Fixing the qubits to be used on the QPU 
layout = [37, 38, 39, 40 ,41, 53, 60, 59, 58, 57, 56, 52]

# Fanne uno per ognuno degli hardware in modo da utilizzare
# sempre il migliore possibile [TODO]

# Puoi anche investigare il plug_in "mapomatic" che vuole
# automatizzare questo step. [TODO]
# https://github.com/qiskit-community/mapomatic

# Trnspiling the circuit to match the backend
pm = generate_preset_pass_manager(optimization_level=2, backend=backend, initial_layout=layout)

# Transpile the circuit to an "Instruction Set Architecture" (ISA) circuit.
# Note: the transpiler automatically adds "ancilla" qubits to make the transpiled
# circuit match the size of the FakeSherbrooke backend.
isa_circ = pm.run(circ)

# Allocating the POVMscheme to use for the expectation value estimation
POVM_scheme = ClassicalShadows(num_qubits=n_spin, measurement_twirl=False, shot_repetitions=1)
POVM_scheme.measurement_circuit.draw("mpl")

# ATTENZIONE:: Da eseguire SOLO per lanciare un nuovo job
# Executing on the choosen backend the POVM
id_file = lancia_POVM(backend, isa_circ, POVM_scheme)