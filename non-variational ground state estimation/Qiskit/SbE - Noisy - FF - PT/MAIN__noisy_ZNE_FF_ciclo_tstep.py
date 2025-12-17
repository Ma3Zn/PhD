from qiskit_aer import AerSimulator
from Ring import *

import numpy as np

# Importante metterla dopo tutti gli altri import
from qiskit.primitives import BackendEstimatorV2 as Estimator

# Ring sites
n_spin = 4

backend_choice = "Fez"

# Setting the range of powers of the Hamiltonian to compute
min_pow = 1
max_pow = 1

# Number of trotter step -- each trotter step is 15-CNOT deep
n_tstep_min = 3
n_tstep_max = 3

# Folding Factors for the ZNE
#folding_factors = np.linspace(1,3,17)
#folding_factors = [1.0, 1.05]
# folding_factors = [1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3]
folding_factors = [1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65]
# folding_factors = [1.75, 1.85, 2, 2.25, 2.5, 2.75]

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

# Selecting the best layout
layout = range(0,n_spin)

# allocazione estimator rumoroso
estimator = Estimator(backend=sim)
estimator.options.default_precision = 0.1

for n_tstep in range(n_tstep_min, n_tstep_max+1):
#   Output avanzamento
    print(f"Inizo simulazione circuito con {n_tstep} Trotter step.")

#   Apertura file output di output
    dir = "./Results/" + backend_choice + "/" + str(n_tstep)

#   Output di DEBUG
    print(f"Output soluzioni al percorso {dir}")

    output_mean = open(dir + "/mean_exp_val.txt", "w")
    output_var  = open(dir + "/var_exp_val.txt", "w")
    output_ff   = open(dir + "/folding_factors.txt", "w")

#   Creation of the annealing circuit
    qc = trotterized_annealing(n_spin, n_tstep)

#   Computing the varios expectation values to different noise strenght levels
    compute_expextation_values_SbE_for_ZNE(qc, H, estimator, min_pow, max_pow,  \
                                           layout, my_backend, folding_factors, \
                                           output_mean, output_var, output_ff)
        
#   Chiudiamo lo straming di output
    output_mean.flush()
    output_mean.close()

    output_var.flush()
    output_var.close()

    output_ff.flush()
    output_ff.close()
