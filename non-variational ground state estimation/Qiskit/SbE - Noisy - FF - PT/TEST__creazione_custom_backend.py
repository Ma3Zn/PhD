from CustomBackend import *
from qiskit_aer import AerSimulator
from qiskit.primitives import BackendEstimatorV2 as Estimator

# Istanzia un backend con 5 qubit in anello
my_backend = MyFakeRingBackend(n_qubits=5,
                                 single_qubit_error=0.002,
                                 single_qubit_duration=50e-9,
                                 two_qubit_error=0.02,
                                 two_qubit_duration=300e-9,
                                 measure_error=0.002)

print("Numero di qubit:", my_backend.num_qubits)
print("Operazioni supportate:", my_backend.target.operation_names)
print(my_backend.layout)

# Possiamo verificare la connettività:
for gate_name in my_backend.target.operation_names:
    print(f"Gate: {gate_name} -> qubit/coppie supportate =",
          list(my_backend.target[gate_name].keys()))

# Se vogliamo solo transpilarci un circuito, lo facciamo così:
from qiskit import QuantumCircuit, transpile

qc = QuantumCircuit(5)
qc.h(0)
qc.rz(1.2,2)
qc.cx(0,1)
qc.cx(1,2)
qc.measure_all()

transpiled_qc = transpile(qc, my_backend)
print("Circuito transpiled:")
print(transpiled_qc)

# creazione di un simulatore rumoroso
noise_model = build_noise_model_from_target(my_backend)
sim = AerSimulator(noise_model=noise_model)

estimator = Estimator(backend=sim)

# Esempio: calcolo <ZZ> su 2 qubit
from qiskit import QuantumCircuit
from qiskit.quantum_info import Pauli

qc = QuantumCircuit(5)
# qc.h(0)
# qc.h(1)
# qc.h(2)
# qc.h(3)
# qc.h(4)
# qc.x(0)
# qc.cx(0, 1)

transpiled_qc = transpile(qc, my_backend)
print("Circuito transpiled:")
print(transpiled_qc)

# Costruiamo un'osservabile ZZ
obs = Pauli("IIIIZ")
print(obs)

# Eseguiamo la stima
job = estimator.run([(transpiled_qc, obs)], precision=1e-3)
result = job.result()[0]
print("Expectation values:", result.data.evs)
