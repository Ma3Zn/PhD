from qiskit_aer import AerSimulator
from Ring import *

import numpy as np

# Importante metterla dopo tutti gli altri import
from qiskit.primitives import BackendEstimatorV2 as Estimator

# Ring sites
n_spin = 8

backend_choice = "Fez"

# Setting the range of powers of the Hamiltonian to compute
min_pow = 4
max_pow = 5

# Number of trotter step -- each trotter step is 15-CNOT deep
n_tstep = 7

# Apertura file output di output
dir = "./Results/Expectation Values - ZNE - FF/" + backend_choice + "/" + str(n_tstep)

output    = open(dir + "/exp_values.txt", "w")
output_ff = open(dir + "/folding_factors.txt", "w")

if backend_choice == "Brisbane":
#   Istanzia un backend con n_spin qubit in anello [Brisbane like]
    my_backend = MyFakeRingBackend(n_qubits=n_spin,
                                    single_qubit_error=2.236e-4,
                                    single_qubit_duration=50e-9,
                                    two_qubit_error=7.519e-3,
                                    two_qubit_duration=660e-9,
                                    measure_error=1.660e-2,
                                    measure_duration=1300)
elif backend_choice == "Fez":
#   Istanzia un backend con n_spin qubit in anello [fez like]
    my_backend = MyFakeRingBackend(n_qubits=n_spin,
                                    single_qubit_error=2.297e-4,
                                    single_qubit_duration=50e-9,
                                    two_qubit_error=4.047e-3,
                                    two_qubit_duration=84e-9,
                                    measure_error=8.179e-3,
                                    measure_duration=1560)
elif backend_choice == "Future":
#   Istanzia un backend con n_spin qubit in anello [future]
    my_backend = MyFakeRingBackend(n_qubits=n_spin,
                                    single_qubit_error=2.297e-5,
                                    single_qubit_duration=50e-9,
                                    two_qubit_error=4.047e-4,
                                    two_qubit_duration=84e-9,
                                    measure_error=1.179e-3,
                                    measure_duration=1560)
else:
# Istanzia un backend con n_spin qubit in anello [IDEALE]
    my_backend = MyFakeRingBackend(n_qubits=n_spin,
                                    single_qubit_error=0,
                                    single_qubit_duration=50e-9,
                                    two_qubit_error=0,
                                    two_qubit_duration=84e-9,
                                    measure_error=0,
                                    measure_duration=1560)

# Creation of the system Hamiltonian as SparsePauliOp
H = build_hamiltonian(n_spin)

# creazione di un simulatore rumoroso
noise_model = build_noise_model_from_target(my_backend)
sim = AerSimulator(noise_model=noise_model)

# allocazione estimator rumoroso
estimator = Estimator(backend=sim)
estimator.options.default_precision = 0.01

# Output avanzamento
print(f"Inizo simulazione circuito con {n_tstep} Trotter step.")

# Creation of the annealing circuit
qc = trotterized_annealing(n_spin, n_tstep)

# Selecting the best layout
layout = range(0,n_spin)

# Folding Factors for the ZNE
# folding_factors = np.linspace(1,3,41)
folding_factors = np.linspace(1,3,17)

# Computing the varios expectation values to different noise strenght levels
compute_expextation_values_SbE_for_ZNE(qc, H, estimator, min_pow, max_pow,  \
                                       layout, my_backend, folding_factors, \
                                       output, output_ff)
        
# Chiudiamo lo straming di output
output.flush()
output.close()