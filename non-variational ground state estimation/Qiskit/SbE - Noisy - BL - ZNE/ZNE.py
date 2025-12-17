import mitiq
from mitiq import zne

from mitiq.zne.scaling import fold_gates_at_random

from qiskit import *
from qiskit.quantum_info import *
from qiskit_ibm_runtime import EstimatorV2 as Estimator
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

from HamOp import *

def fold_circuit(circuit: QuantumCircuit,   \
                 scale_factor: list[int],   \
                 layout: list[int],         \
                 backend):
    
#   Ritorna il circuito ampliato per un dato fattore di scala.
#   Se scale_factor == 1, ritorna il circuito originale.
#   Per scale_factor > 1, fold_gates_at_random aggiunge coppie di porte
#   che si annullano idealmente, ma che amplificano il rumore.
    
    if scale_factor == 1:
        qc = circuit
    else:
        qc = fold_gates_at_random(circuit, scale_factor)

#   Eseguiamo un transpiling del circuito per matchare il backend
    pm = generate_preset_pass_manager(optimization_level=0, backend=backend, initial_layout=layout)

#   Transpile the circuit to an "Instruction Set Architecture" (ISA) circuit.
#   Note: the transpiler automatically adds "ancilla" qubits to make the transpiled
#   circuit match the size of the backend.
    isa_circ = pm.run(qc)

    return isa_circ

def stima_aspettazione_ZNE(circuit: QuantumCircuit,     \
                            oss: SparsePauliOp,         \
                            estimator: Estimator,       \
                            fold_factors: list[int],    \
                            layout: list[int],          \
                            backend,                    \
                            out) -> float:
#   Funzione che dato un circuito, un'osservabile ed un estimator, stima il 
#   valore di aspettazione di questa eseguendo un'estrapolazione dell'errore

#   Creiamo una lista vuota che conterra' i valori dell'osservabile per i vari 
#   livelli di folding
    exp_values = []

#   Cicliamo sui vari fattori di scala per generare gli opportuni circuiti
#   ampliati
    for ff in fold_factors:

#       Output di avanzamento
        print("Fold Factor: %.1f" % ff)

#       Creiamo il circuito ampliato
        isa_folded_circuit = fold_circuit(circuit, ff, layout, backend)

#       Applichiamo il layout opportuno all'osservabile
        isa_oss = oss.apply_layout(isa_folded_circuit.layout)

#       Creazione ed esecuzione del job
        job = estimator.run([(isa_folded_circuit, isa_oss)])
        pub_result = job.result()[0]

#       Salvataggio valore d'aspettazione
        exp_values.append(pub_result.data.evs)

#       Stampiamo su file tutti i valori dell'osservabile attuale per poi poter 
#       verificare se la base scelta per eseguire il fit è la miglior scelta o meno
        out.write("%.16e " % pub_result.data.evs)
        out.flush()

#   end ff

#   Aggiorniamo lo stream su file
    out.write("\n")

#   ATTENZIONE::    Da decommentare se si vuole tornare già una stima 
#                   del valore di aspettazione "Zero Noised"

# #   Eseguiamo un fit polinomiale di grado 1 [lineare]
#     fit_coeff = np.polyfit(fold_factors, exp_values, 1)

# #   Ritorniamo il valore d'aspettazione dell'osservabile estrapolato a "zero rumore"
#     return fit_coeff[1]

# #   Estrapolazione esponenziale con mitiq
#     return zne.ExpFactory.extrapolate(fold_factors, exp_values, asymptote=0.5)