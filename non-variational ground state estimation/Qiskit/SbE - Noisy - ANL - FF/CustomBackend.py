from qiskit.providers.backend import BackendV2
from qiskit.providers import Options
from qiskit.transpiler import Target, InstructionProperties
from qiskit.circuit.library.standard_gates import U3Gate, CXGate
from qiskit.circuit import Instruction

from qiskit_aer.noise import NoiseModel, errors

class MyFakeRingBackend(BackendV2):
    def __init__(
        self,
        n_qubits: int,
        # Errori parametrici
        single_qubit_error: float = 0.001,
        two_qubit_error: float = 0.01,
        measure_error: float = 0.02,
        # Durate parametriche (in secondi)
        single_qubit_duration: float = 50e-9,
        two_qubit_duration: float = 300e-9,
        measure_duration: float = 1e-6,
        # Risoluzione temporale (per scheduling)
        dt: float = 1e-9,
    ):
        """
        Backend V2 fittizio con topologia ad anello e parametri d'errore/durata parametrici.
        
        Args:
            n_qubits: Numero di qubit.
            single_qubit_error: Errore (approssimato) per gate 1-qubit (u3).
            two_qubit_error: Errore (approssimato) per gate 2-qubit (cx).
            measure_error: Errore (approssimato) di misura.
            single_qubit_duration: Durata gate 1-qubit (in secondi).
            two_qubit_duration: Durata gate 2-qubit (in secondi).
            measure_duration: Durata misura (in secondi).
            dt: Valore di time step nel Target (usato dal transpiler in alcune pass).
        """
        super().__init__(name=f"MyFakeRingBackend_{n_qubits}")
        self._n_qubits = n_qubits
        self._options = Options()
        self._options.shots = 1024

        # 1) Creiamo un Target, che includa dt
        self._target = Target(num_qubits=n_qubits, dt=dt)

        # 2) Aggiungiamo gate single-qubit (u3) su tutti i qubit
        u3_map = {}
        for q in range(n_qubits):
            u3_map[(q,)] = InstructionProperties(
                duration=single_qubit_duration,
                error=single_qubit_error
            )
        self._target.add_instruction(U3Gate(0, 0, 0), u3_map)

        # 3) Aggiungiamo gate two-qubit (cx) ad anello, bidirezionale
        cx_map = {}
        for i in range(n_qubits):
            j = (i + 1) % n_qubits
            cx_map[(i, j)] = InstructionProperties(
                duration=two_qubit_duration,
                error=two_qubit_error
            )
            cx_map[(j, i)] = InstructionProperties(
                duration=two_qubit_duration,
                error=two_qubit_error
            )
        self._target.add_instruction(CXGate(), cx_map)

        # 4) Aggiungiamo un'istruzione "measure" (fittizia) su ogni qubit
        measure_instr = Instruction(
            name="measure",
            num_qubits=1,
            num_clbits=1,
            params=[]
        )
        measure_map = {}
        for q in range(n_qubits):
            measure_map[(q,)] = InstructionProperties(
                duration=measure_duration,
                error=measure_error
            )
        self._target.add_instruction(measure_instr, measure_map)

        # Parametri di "rumore" memorizzati (non strettamente necessari, ma comodi).
        self._single_qubit_error = single_qubit_error
        self._two_qubit_error = two_qubit_error
        self._measure_error = measure_error

    @property
    def num_qubits(self) -> int:
        return self._n_qubits

    @property
    def target(self) -> Target:
        return self._target

    def run(self, circuits, **run_options):
        """
        Potresti implementare un'esecuzione reale (delegata ad AerSimulator) qui,
        ma in questo esempio solleviamo un'eccezione.
        """
        raise NotImplementedError("MyFakeRingBackend non esegue direttamente i circuiti.")
    
    def max_circuits(self, **run_options) -> int:
        return 999999

    @classmethod
    def _default_options(cls):
        return Options(shots=1024)

def build_noise_model_from_target(backend: MyFakeRingBackend) -> NoiseModel:
    noise_model = NoiseModel()
    for gate_name in backend.target.operation_names:
        inst_map = backend.target[gate_name]
        for qubits, props in inst_map.items():
            err = props.error
            if err is not None and err > 0:
                if gate_name == "measure":
                    # Usiamo un readout error 1-qubit (assumendo un solo qubit misurato)
                    # Struttura: [[1 - e, e],[e, 1 - e]]
                    readout_error = errors.ReadoutError([[1-err, err],
                                                         [err,   1-err]])
                    noise_model.add_readout_error(readout_error, qubits)
                else:
                    # Gate quantistici (u3, cx) -> canale depolarizzante
                    num_q = len(qubits)
                    depol = errors.depolarizing_error(err, num_q)
                    noise_model.add_quantum_error(depol, gate_name, qubits)
    return noise_model