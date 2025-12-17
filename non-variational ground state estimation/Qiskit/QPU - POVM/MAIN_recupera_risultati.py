from qiskit.primitives import Estimator

import scipy.linalg as ln

from povm_toolbox.library import *

from Ring import *

from POVM import recupera_POVM

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
n_tstep = 5

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

#ATTENZIONE:: Da eseguire SOLO per recuperare i risultati di un job terminato
# Retriving the job results
id_file = '10000_5stp_dd'
POVM_results = recupera_POVM(service, id_file)

# Setting the maximum power of the Hamiltonian to use
max_pow = 1

# Mappare le misure delle osservabili al layout del cricuito attuale
# non credo sia da fare. Credo che sia gestito in automatico, ma va
# sicuramente indagato in merito [TODO]

# Output avanzamento
print("Risultati recuperati -- Inizio stima parametri")
print(id_file)

# Retriving an energy estimation throught digitalized annealing + subspace expansion
estimate = energy_estimate_SbE(isa_circ, H, POVM_results, max_pow)

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