from qiskit_aer import AerSimulator

from Ring import *

# Importante metterla dopo tutti gli altri import
from qiskit_ibm_runtime import QiskitRuntimeService, EstimatorV2
from qiskit_ibm_runtime.noise_learner import NoiseLearner
from qiskit_ibm_runtime.options import (
    NoiseLearnerOptions,
    ResilienceOptionsV2,
    EstimatorOptions,
)

# Ring sites
n_spin = 8

# Creazione layout fittizio backend
layout = range(0, n_spin)

# Istanzia un backend con n_spin qubit in anello [Brisbane like]
my_backend = MyFakeRingBackend(n_qubits=n_spin,
                                 single_qubit_error=2.236e-4,
                                 single_qubit_duration=50e-9,
                                 two_qubit_error=7.519e-3,
                                 two_qubit_duration=660e-9,
                                 measure_error=1.660e-2,
                                 measure_duration=1300)

# # Istanzia un backend con n_spin qubit in anello [fez like]
# my_backend = MyFakeRingBackend(n_qubits=n_spin,
#                                  single_qubit_error=2.297e-4,
#                                  single_qubit_duration=50e-9,
#                                  two_qubit_error=4.047e-3,
#                                  two_qubit_duration=84e-9,
#                                  measure_error=8.179e-3,
#                                  measure_duration=1560)

# # Istanzia un backend con n_spin qubit in anello [less far future]
# my_backend = MyFakeRingBackend(n_qubits=n_spin,
#                                  single_qubit_error=2.297e-5,
#                                  single_qubit_duration=50e-9,
#                                  two_qubit_error=4.047e-4,
#                                  two_qubit_duration=84e-9,
#                                  measure_error=1.179e-3,
#                                  measure_duration=1560)

# # Istanzia un backend con n_spin qubit in anello [far future]
# my_backend = MyFakeRingBackend(n_qubits=n_spin,
#                                  single_qubit_error=2.297e-6,
#                                  single_qubit_duration=50e-9,
#                                  two_qubit_error=4.047e-5,
#                                  two_qubit_duration=84e-9,
#                                  measure_error=1.179e-4,
#                                  measure_duration=1560)

# # Istanzia un backend con n_spin qubit in anello [IDEALE]
# my_backend = MyFakeRingBackend(n_qubits=n_spin,
#                                  single_qubit_error=0,
#                                  single_qubit_duration=50e-9,
#                                  two_qubit_error=0,
#                                  two_qubit_duration=84e-9,
#                                  measure_error=0,
#                                  measure_duration=1560)

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

# creazione di un simulatore rumoroso
noise_model = build_noise_model_from_target(my_backend)
sim = AerSimulator(noise_model=noise_model)

# Setting the maximum power of the Hamiltonian to use
max_pow = 1

# Apertura file output .txt
output = open("./Results/err.txt", "w")

# Allocazione Pass Menager - Transpiler
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
pm = generate_preset_pass_manager(optimization_level=2, backend=my_backend)

# Cicliamo sui vari trotter step da eseguire
for n_tstep in range(0,n_tstep_max + 1):
#   Output avanzamento
    print(f"Inizo simulazione circuito con {n_tstep} Trotter step.")

#   Creation of the annealing circuit
    qc = trotterized_annealing(n_spin, n_tstep)

#   Transpile the circuit to an "Instruction Set Architecture" (ISA) circuit.
#   Note: the transpiler automatically adds "ancilla" qubits to make the transpiled
#   circuit match the size of the FakeSherbrooke backend.
    isa_circ = pm.run(qc)

# #   Instantiate a noise learner options object
#     learner_options = NoiseLearnerOptions(
#         max_layers_to_learn=15, num_randomizations=32, twirling_strategy="all"
#     )
 
# #   Instantiate a NoiseLearner object and execute the noise learning program
#     learner = NoiseLearner(mode=sim, options=learner_options)
#     job = learner.run([isa_circ])
#     noise_model = job.result()

#   allocazione estimator rumoroso
    estimator = Estimator(mode=sim)
    estimator.options.default_precision = 0.01

#   abilitiamo la mitigazione degli errori tramite PEA - Pauli Twirling - TREX
    # estimator.options.resilience.layer_noise_model = noise_model

    estimator.options.twirling.enable_gates = True
    estimator.options.twirling.num_randomizations = 32
    estimator.options.twirling.shots_per_randomization = 100

    estimator.options.resilience.measure_mitigation = True
    estimator.options.resilience.measure_noise_learning.num_randomizations = 32
    estimator.options.resilience.measure_noise_learning.shots_per_randomization = 100

    estimator.options.resilience.zne_mitigation = True
    estimator.options.resilience.zne.amplifier = "pea"

#   ATTENZIONE:: con questa logica sprechiamo enorme tempo nel ricalcolare per ogni 
#   circuito tutte le potenze dell'Hamiltoniana, che sono INDIPENDENTI dal circuito.
#   In una versione futura di queste simulazioni è tassativo invertire l'ordine delle
#   operazioni, il guadagno sarà notevole.

#   Retriving an energy estimation throught digitalized annealing + subspace expansion
    try:
        estimate = energy_estimate_SbE(isa_circ, H, estimator, max_pow, \
                                       layout, my_backend, output)
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