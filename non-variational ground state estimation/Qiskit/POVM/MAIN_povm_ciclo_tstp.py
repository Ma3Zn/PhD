from qiskit.primitives import Estimator

import scipy.linalg as ln

from povm_toolbox.library import *

from Ring import *

# Ring sites
n_spin = 4

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

# Setting the maximum power of the Hamiltonian to use
max_pow = 1

# Apertura file output .txt
output = open("./Results/err.txt", "w")

# Cicliamo sui vari trotter step da eseguire
for n_tstep in range(0,n_tstep_max + 1):
#   Output avanzamento
    print(f"Inizo simulazione circuito con {n_tstep} Trotter step.")

#   Creation of the annealing circuit
    circ = trotterized_annealing(n_spin, n_tstep)

#   Allocating the POVMscheme to use for the expectation value estimation
    POVM_scheme = ClassicalShadows(num_qubits=n_spin, measurement_twirl=False, shot_repetitions=1)
    POVM_scheme.measurement_circuit.draw("mpl")

#   Retriving an energy estimation throught digitalized annealing + subspace expansion
    try:
        estimate = energy_estimate_SbE(circ, H, POVM_scheme, max_pow)
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