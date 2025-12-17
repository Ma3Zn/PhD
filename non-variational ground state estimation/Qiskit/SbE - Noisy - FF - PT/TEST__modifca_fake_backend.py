from qiskit_ibm_runtime.fake_provider import FakeTorino
from qiskit.providers.backend import BackendV2
from qiskit.transpiler import Target, InstructionProperties
from qiskit_aer.noise import NoiseModel

# -----------------------------------------------------------------------------
# STEP 1: Estrazione del fake backend e dei dati utili
# -----------------------------------------------------------------------------
fake_backend = FakeTorino()
noise_model_full = NoiseModel.from_backend(fake_backend)
try:
    gate_lengths = fake_backend.configuration().gate_lengths
except Exception:
    gate_lengths = {}

# -----------------------------------------------------------------------------
# STEP 2: Creazione di un nuovo backend con N qubit e topologia "ring"
# -----------------------------------------------------------------------------
N = 5  # Numero di qubit del nuovo backend
num_qubits = N

# Costruiamo la coupling map "ring" bidirezionale:
ring_coupling_map = []
for i in range(num_qubits):
    j = (i + 1) % num_qubits
    ring_coupling_map.append((i, j))
    ring_coupling_map.append((j, i))

# -----------------------------------------------------------------------------
# STEP 3: Creazione del nuovo Target con le istruzioni e le proprietà
# -----------------------------------------------------------------------------
new_target = Target(num_qubits=num_qubits)

def safe_add_instruction(target, instr, mapping):
    """Aggiunge mapping all'istruzione; se l'istruzione è già presente, aggiorna la mappa."""
    if instr in target.instructions:
        target.instructions[instr].update(mapping)
    else:
        target.add_instruction(instr, mapping)

def get_error_for_gate(gate_name, qubits):
    """
    Recupera l'errore per il gate (gate_name) sui qubit specificati.
    
    Cerca la chiave (gate_name, tuple(qubits)) nel NoiseModel.
    Se non viene trovato e l'operazione è a due qubit, restituisce il primo errore
    trovato per quel gate (tra quelli disponibili).
    """
    key = (gate_name, tuple(qubits))
    errors = noise_model_full._local_quantum_errors.get(key, None)
    if errors is not None and len(errors) > 0:
        return errors[0]
    if len(qubits) == 2:
        # Se non troviamo l'errore specifico, iteriamo sulle chiavi disponibili a 2 qubit
        for k, err_list in noise_model_full._local_quantum_errors.items():
            if isinstance(k, tuple) and len(k) == 2:
                g_name, qtuple = k
                if g_name == gate_name and len(qtuple) == 2:
                    return err_list[0]
    return None

def get_gate_time(gate_name, qubits):
    """Recupera il tempo di gate dalla configurazione, oppure usa un valore di default."""
    default = 50 if len(qubits) == 1 else 200
    for key, t in gate_lengths.items():
        if key[0] == gate_name:
            return t
    return default

# Iteriamo sulle istruzioni del fake backend.
# Nota: fake_backend.target.instructions è una lista di tuple (gate_obj, mapping_dict)
for gate_obj, instr_props in fake_backend.target.instructions:
    if gate_obj.num_qubits == 1:
        # Per ogni qubit del nuovo backend aggiungiamo l'operazione a singolo qubit.
        for qubit in range(num_qubits):
            error = get_error_for_gate(gate_obj.name, (qubit,))
            gate_time = get_gate_time(gate_obj.name, (qubit,))
            ip = InstructionProperties()
            ip.error = error
            ip.duration = gate_time
            safe_add_instruction(new_target, gate_obj, {(qubit,): ip})
    elif gate_obj.num_qubits == 2:
        # Per operazioni a due qubit, definiamo le coppie uniche dalla nostra topologia "ring".
        unique_edges = set()
        for edge in ring_coupling_map:
            sorted_edge = tuple(sorted(edge))
            unique_edges.add(sorted_edge)
        for edge in unique_edges:
            error = get_error_for_gate(gate_obj.name, edge)
            if error is None:
                # Se non è disponibile un errore per questa coppia, proviamo a usare uno preso dalle proprietà originali
                for prop in instr_props.values():
                    error = prop.get('error', None)
                    if error is not None:
                        break
            gate_time = get_gate_time(gate_obj.name, edge)
            ip = InstructionProperties()
            ip.error = error
            ip.duration = gate_time
            safe_add_instruction(new_target, gate_obj, {edge: ip})

# -----------------------------------------------------------------------------
# STEP 4: Definizione del nuovo backend basato su BackendV2
# -----------------------------------------------------------------------------
class FakeRingBackend(BackendV2):
    def __init__(self, num_qubits, target, coupling_map):
        self._num_qubits = num_qubits
        self._target = target
        self._coupling_map = coupling_map

    @property
    def backend_name(self):
        return "FakeRingBackend"

    @property
    def num_qubits(self):
        return self._num_qubits

    @property
    def target(self):
        return self._target

    @property
    def coupling_map(self):
        return self._coupling_map

    def run(self, circuits, **run_options):
        raise NotImplementedError("Questo fake backend è solo per scopi di configurazione.")

fake_ring_backend = FakeRingBackend(num_qubits, new_target, ring_coupling_map)

# -----------------------------------------------------------------------------
# STEP 5: Verifica
# -----------------------------------------------------------------------------
print("Backend Name:", fake_ring_backend.backend_name)
print("Numero di Qubits:", fake_ring_backend.num_qubits)
print("Coupling Map (ring):", fake_ring_backend.coupling_map)
print("Target instructions:")
for gate_obj, mapping in fake_ring_backend.target.instructions.items():
    print(f"Gate: {gate_obj.name}, mapping: {mapping}")
