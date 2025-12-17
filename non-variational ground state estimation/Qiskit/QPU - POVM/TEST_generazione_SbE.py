from qiskit.primitives import Estimator

import scipy.linalg as ln

from povm_toolbox.library import *

from Ring import *

from qiskit_ibm_runtime import QiskitRuntimeService

service = QiskitRuntimeService(
    channel='ibm_quantum',
    instance='ibm-q/open/main',
    token='feb77ce69166caff54852530c1634b44e307ee4fd5af82cbea78c06774b07ed575a010878366428456e3333366c38aed28bc7a6dd61279eca7febcebc353f3b3'
)

# Ring sites
n_spin = 12

# Creation of the system Hamiltonian as SparsePauliOp
H = build_hamiltonian(n_spin)

# Number of trotter step -- each trotter step is 15-CNOT deep
n_tstep = 7

# Creation of the annealing circuit
qc = trotterized_annealing(n_spin, n_tstep)

backend = service.backend("ibm_brisbane")

print(type(backend))
print(backend)

# Fixing the qubits to be used on the QPU 
layout = [37, 38, 39, 40 ,41, 53, 60, 59, 58, 57, 56, 52]

from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
pm = generate_preset_pass_manager(optimization_level=2, backend=backend, initial_layout=layout)

# Transpile the circuit to an "Instruction Set Architecture" (ISA) circuit.
# Note: the transpiler automatically adds "ancilla" qubits to make the transpiled
# circuit match the size of the FakeSherbrooke backend.
isa_circ = pm.run(qc)

# circ.draw()

from POVM import recupera_POVM
#ATTENZION:: Da eseguire SOLO per recuperare i risultati di un job terminato

id_file = '20000_2024-10-28 13:47:26.509748'

# Retriving the job results
POVM_results = recupera_POVM(service, id_file)

# Setting the maximum power of the Hamiltonian to use
max_pow = 1

# Mappare le misure delle osservabili al layout del cricuito attuale [TODO]
# Retriving an energy estimation throught digitalized annealing + subspace expansion
estimate = energy_estimate_SbE(isa_circ, H, POVM_results, max_pow)