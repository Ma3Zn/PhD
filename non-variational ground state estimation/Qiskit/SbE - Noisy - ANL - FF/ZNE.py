import mitiq
from mitiq import zne

from mitiq.zne.scaling import fold_gates_at_random

from qiskit import *
from qiskit.quantum_info import *
from qiskit_ibm_runtime import EstimatorV2 as Estimator
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

from HamOp import *

def count_gates(circuit: QuantumCircuit):
#   Funzione che calcola il numero di gate nel circuito passato
    gates = dict(circuit.count_ops())
    n_gates = 0
    for x in gates:
        n_gates = n_gates + gates[x]

    print(gates)

    return n_gates

def fold_circuit(circuit: QuantumCircuit,   \
                 scale_factor: list[int],   \
                 layout: list[int],         \
                 backend):
    
#   Ritorna il circuito ampliato per un dato fattore di scala.
#   Se scale_factor == 1, ritorna il circuito originale.
#   Per scale_factor > 1, fold_gates_at_random aggiunge coppie di porte
#   che si annullano idealmente, ma che amplificano il rumore.
    
#   Calcoliamo il numero di gate iniziali
    n_gate_init = count_gates(circuit)

    if scale_factor == 1:
        qc = circuit
    else:
        qc = fold_gates_at_random(circuit, scale_factor, fidelities={"single": 1.0, "CNOT": 0.9960})

#   Eseguiamo un transpiling del circuito per matchare il backend
    pm = generate_preset_pass_manager(optimization_level=0, backend=backend, initial_layout=layout)

#   Transpile the circuit to an "Instruction Set Architecture" (ISA) circuit.
#   Note: the transpiler automatically adds "ancilla" qubits to make the transpiled
#   circuit match the size of the backend.
    isa_circ = pm.run(qc)

#   Contiamo il numero di gate totali
    n_gate_fold = count_gates(isa_circ)

#   Calcoliamo l'effettivo folding_factor
    ff_effettivo = n_gate_fold / n_gate_init

    return isa_circ, ff_effettivo

def stima_aspettazione_ZNE(circuit: QuantumCircuit,     \
                            oss: SparsePauliOp,         \
                            estimator: Estimator,       \
                            fold_factors: list[int],    \
                            layout: list[int],          \
                            backend,                    \
                            out, out_ff) -> float:
#   Funzione che dato un circuito, un'osservabile ed un estimator, stima il 
#   valore di aspettazione di questa eseguendo un'estrapolazione dell'errore

#   Inizializziamo il vettore che conterra' i folding factors effettivi
    fold_factors_effettivi = []

#   Inizializziamo la lista di osservabili da misurare
    isa_oss = []

#   Applichiamo il layout opportuno alle varie osservabili
    for o in oss:
        isa_oss.append(o.apply_layout(layout))

#   Cicliamo sui vari fattori di scala per generare gli opportuni circuiti
#   ampliati
    for ff in fold_factors:

#       Creiamo il circuito ampliato
        isa_folded_circuit, ff_effettivo = fold_circuit(circuit, ff, layout, backend)

#       Output di avanzamento
        print("Fold Factor: %.3f" % ff_effettivo)

#       Scriviamo su file il fold_factor effettivo
        out_ff.write("%.16f " % ff_effettivo)
        out_ff.flush()

#       Creazione ed esecuzione del job
        job = estimator.run([(isa_folded_circuit, isa_oss)])
        pub_result = job.result()[0]

#       Stampiamo su file tutti i valori dell'osservabile attuale per poi poter 
#       verificare se la base scelta per eseguire il fit è la miglior scelta o meno
        for exp_val in pub_result.data.evs:
            out.write("%.16e " % exp_val)

#       Aggiorniamo lo stream su file
        out.write("\n")
        out.flush()

#   end ff
