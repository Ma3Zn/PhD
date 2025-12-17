from qiskit_ibm_runtime.fake_provider import FakeBrisbane
from qiskit_ibm_runtime import EstimatorV2 as Estimator

import mapomatic as mm

from Ring import *


fake_backend = FakeBrisbane()

# Ring sites
n_spin = 12

# Creation of the system Hamiltonian as SparsePauliOp
H = build_hamiltonian(n_spin)

# Number of trotter step -- each trotter step is 15-CNOT deep
n_tstep = 4

# Allocazione Estimator da usare nella stima dei parametri della SbE
# estimator = Estimator() # EXACT Estimator
estimator = Estimator(mode=fake_backend)

# Setting the maximum power of the Hamiltonian to use for the Krylov space
max_pow = 3

# Apertura file output .txt
output = open("./Results/exp_values.txt", "w")

# Output avanzamento
print(f"Inizo simulazione circuito con {n_tstep} Trotter step.")

# Creation of the annealing circuit
qc = trotterized_annealing(n_spin, n_tstep)

# Transpiling the circuit for the fake backend
trans_qc = transpile(qc, fake_backend, optimization_level=3)

# Reducing the size of the circuit
small_qc = mm.deflate_circuit(trans_qc)

# Finding the matchin subgraph for the fake backend
layouts = mm.matching_layouts(small_qc, fake_backend)

# Evaluating the cost of each layout
scores = mm.evaluate_layouts(small_qc, layouts, fake_backend)

# Selecting the best layout
layout = scores[0][0]

# Folding Factors for the ZNE
folding_factors = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]

# Computing the varios expectation values to different noise strenght levels
compute_expextation_values_SbE_for_ZNE(qc, H, estimator, max_pow, \
                                       layout, fake_backend,      \
                                       folding_factors, output)
        
# Chiudiamo lo straming di output
output.flush()
output.close()