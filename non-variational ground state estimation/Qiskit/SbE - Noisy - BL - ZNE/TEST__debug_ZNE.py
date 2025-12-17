from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel
# from qiskit.primitives import BackendEstimator as Estimator

from qiskit_ibm_runtime import EstimatorV2 as Estimator

import scipy.linalg as ln

from Ring import *

from qiskit_ibm_runtime.fake_provider import FakeBrisbane

fake_backend = FakeBrisbane()

# Fixing the qubits to be used on the QPU 
layout = [37, 38, 39, 40 ,41, 53, 60, 59, 58, 57, 56, 52]

# Ring sites
n_spin = 12

# Creation of the system Hamiltonian as SparsePauliOp
H = build_hamiltonian(n_spin)

# Number of trotter step -- each trotter step is 15-CNOT deep
n_tstep = 4

# Creation of the annealing circuit
qc = trotterized_annealing(n_spin, n_tstep)

# Setting the maximum power of the Hamiltonian to use
max_pow = 1

# Creazione del simulatore rumoroso
# Utilizziamo un noise model estratto da un backend fake (FakeVigo)
# noise_model = NoiseModel.from_backend(fake_backend)
# simulator = AerSimulator(noise_model=noise_model)

# # Allocating the estiamtor to use for the Subpspace Expansion -- EXACT Estimator
# estimator = Estimator()

# # Da usare con AerSimulator
# estimator = Estimator(simulator)

estimator = Estimator(mode=fake_backend)

# Apriamo il file di output per la stampa dei valori di aspettazione calcolati
output = open("./Results/exp_values.txt", "w")

# Retriving an energy estimation throught digitalized annealing + subspace expansion
print("Inizio SbE con Zero noise Extrapolation")
estimate = energy_estimate_SbE(qc, H, estimator, max_pow, layout, fake_backend, output)

# Chiudiamo il file di output per la stampa dei valori di aspettazione calcolati
output.close()

# Computing a numerical value of the system Hamiltonian [TBD]
H_num = evaluate_hamiltonian(H)

# Numerical diagonalization of the Hamiltonian [TBD]
e, V = ln.eigh(np.real(H_num), subset_by_index=[0, 0])

# Recuperiamo l'indice del groung state
id_gs = np.argmin(e)

# Recuperiamo il ground state e la relativa energia
gs     = V[:,id_gs]
energy = e[id_gs]