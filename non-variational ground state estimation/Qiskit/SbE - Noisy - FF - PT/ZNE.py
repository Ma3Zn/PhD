import mitiq
from mitiq import zne

from mitiq.zne.scaling import fold_gates_at_random

from qiskit import *
from qiskit.quantum_info import *
from qiskit_ibm_runtime import EstimatorV2 as Estimator
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

from mitiq.pt import generate_pauli_twirl_variants

from HamOp import *

import numpy as np

def calcola_media_e_varianza(exp_val: np.matrix, n_sample: int):
#   Funzione che passata una matrice avente per ogni colonna i vari di 
#   aspettazione di una determinata osservabile calcolati per diversi twirled 
#   circuit ne calcola la media e la varianza campionaria di questi

#   Calcoliamo la media dei valori di aspettazione (quindi la media a colonna
#   fissata e su tutte le righe)
    media = np.mean(exp_val, axis=0)

#   Calcoliamo la varianza campionaria 
    var = np.sum(np.square(exp_val-media), axis=0) / (n_sample - 1)

    return media, var

def count_gates(circuit: QuantumCircuit):
#   Funzione che calcola il numero di gate nel circuito passato
    gates = dict(circuit.count_ops())
    n_gates = 0
    for x in gates:
        if x == 'cx':
            n_gates = n_gates + gates[x]

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
                            out_mean, out_var, out_ff) -> float:
#   Funzione che dato un circuito, un'osservabile ed un estimator, stima il 
#   valore di aspettazione di questa eseguendo un'estrapolazione dell'errore

#   Inizializziamo la lista di osservabili da misurare
    isa_oss = []

#   Applichiamo il layout opportuno alle varie osservabili
    for o in oss:
        isa_oss.append(o.apply_layout(layout))

#   Cicliamo sui vari fattori di scala per generare gli opportuni circuiti
#   ampliati
    for ff in fold_factors:

#       Inizializziamo la lista che conterra' i valori di aspettazione calcolati
#       per le varie versioni twirled dei circuiti
        exp_val = []

#       Creiamo il circuito ampliato
        isa_folded_circuit, ff_effettivo = fold_circuit(circuit, ff, layout, backend)

#       Output di avanzamento
        print("Fold Factor: %.3f" % ff_effettivo)

#       Scriviamo su file il fold_factor effettivo
        out_ff.write("%.16f " % ff_effettivo)
        out_ff.flush()

# ###############################################################################
# #                            TEST SENZA TWIRLING                              #
# ###############################################################################

#         job = estimator.run([(isa_folded_circuit, isa_oss)])
#         pub_result = job.result()[0]

#         mean_exp_val = np.matrix(pub_result.data.evs)
#         var_exp_val  = np.matrix(0)

# ###############################################################################
# #                            TEST SENZA TWIRLING                              #
# ###############################################################################

#       Definiamo il numero di circuiti Twirled da utilizzare
        n_twirl = 30

#       Generiamo i circuiti Twirled
        twirled_circuits = generate_pauli_twirl_variants(isa_folded_circuit, num_circuits=n_twirl)

#       Cicliamo sui circuiti Twirled
        for circ in twirled_circuits:
#           Creazione ed esecuzione del job
            job = estimator.run([(circ, isa_oss)])
            pub_result = job.result()[0]

#           Recuperiamo i valori di aspettazione delle varie osservabili
            exp_val.append(pub_result.data.evs)

#       Calcoliamo media e varianza campionaria dei vari valori di aspettazione
        mean_exp_val, var_exp_val = calcola_media_e_varianza(np.matrix(exp_val), n_twirl)

#       Scriviamo su file le medie
        for i in range(0, len(isa_oss)):
            out_mean.write("%.16e " % mean_exp_val[0,i])

#       Aggiorniamo lo stream su file
        out_mean.write("\n")
        out_mean.flush()

#       Scriviamo su file le varianze
        for i in range(0, len(isa_oss)) :
            out_var.write("%.16e " % var_exp_val[0,i])

#       Aggiorniamo lo stream su file
        out_var.write("\n")
        out_var.flush()

#   end ff
