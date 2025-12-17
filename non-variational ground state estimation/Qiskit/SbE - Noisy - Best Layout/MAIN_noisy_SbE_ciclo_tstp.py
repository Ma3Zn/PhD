from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_ibm_runtime import EstimatorV2 as Estimator
from qiskit_ibm_runtime import SamplerV2

import scipy.linalg as ln

import mapomatic as mm

from Ring import *

from qiskit_ibm_runtime.fake_provider import FakeBrisbane, FakeTorino

fake_backend = FakeTorino()

# Ring sites
n_spin = 12

# Creation of the system Hamiltonian as SparsePauliOp
H = build_hamiltonian(n_spin)

# Computing a numerical value of the system Hamiltonian
H_num = evaluate_hamiltonian(H)

# Numerical diagonalization of the Hamiltonian
e, V = ln.eigh(np.real(H_num), subset_by_index=[0, 0])

# Recuperiamo l'indice del groung state
id_gs = np.argmin(e)

# Recuperiamo il ground state e la relativa energia
gs     = V[:,id_gs]
energy = e[id_gs]

# Number of trotter step -- each trotter step is 15-CNOT deep
n_tstep_max = 7

# Allocazione variabile contenente gli errori relativi per ogni circuito
err = 0.0

# Allocazione Estimator da usare nella stima dei parametri della SbE
# estimator = Estimator() # EXACT Estimator
estimator = Estimator(mode=fake_backend)

# Setting the maximum power of the Hamiltonian to use
max_pow = 1

# Apertura file output .txt
output = open("./Results/err.txt", "w")

# Cicliamo sui vari trotter step da eseguire
for n_tstep in range(0,n_tstep_max + 1):
#   Output avanzamento
    print(f"Inizo simulazione circuito con {n_tstep} Trotter step.")

#   Creation of the annealing circuit
    qc = trotterized_annealing(n_spin, n_tstep)

#   Transpiling the circuit for the fake backend
    trans_qc = transpile(qc, fake_backend, optimization_level=3)

#   Reducing the size of the circuit
    small_qc = mm.deflate_circuit(trans_qc)

#   Finding the matchin subgraph for the fake backend
    layouts = mm.matching_layouts(small_qc, fake_backend)

#   Evaluating the cost of each layout
    scores = mm.evaluate_layouts(small_qc, layouts, fake_backend)

#   Creation of the standard pass menager for the transpiling pipeline with the 
#   best possible initial layout
    pm = generate_preset_pass_manager(optimization_level=3, \
                                      backend=fake_backend, \
                                      initial_layout=scores[0][0])

#   Transpile the circuit to an "Instruction Set Architecture" (ISA) circuit.
#   Note: the transpiler automatically adds "ancilla" qubits to make the
#   transpiled circuit match the size of the FakeSherbrooke backend.
    isa_circ = pm.run(small_qc)

#   Retriving an energy estimation throught digitalized annealing + subspace expansion
    try:
        estimate = energy_estimate_SbE(isa_circ, H, estimator, max_pow)
        err = abs(energy - estimate[0]) / abs(energy)
    except:
        err = -1.0

#   Output finale degli errori su file
    output.write("%.16e " % err)
    output.flush()

#   Aggiornamento avanzamento
    print([energy, estimate, err])

# Chiudiamo lo straming di output
output.close()