from qiskit import *
from qiskit.providers import *
from qiskit.primitives import *
from qiskit.quantum_info import *

from povm_toolbox.library import *
from povm_toolbox.sampler import *
from povm_toolbox.post_processor import *

from povm_toolbox.sampler import POVMSamplerJob

from qiskit_ibm_runtime import SamplerV2 as RuntimeSampler
from qiskit_ibm_runtime import QiskitRuntimeService

from qiskit_ibm_runtime.ibm_backend import IBMBackend

from datetime import datetime

def alloca_POVM_post_processor(result: POVMPubResult) -> POVMPostProcessor:
#   Funzione che alloca il POVMPostProcessor necessario per eseguire la stima
#   dei valori di aspettazione delle osservabili una volta che abbiamo eseguito
#   le POVM. In particolare setiamo il dual frame del post-processor sulla
#   strategia delle frequenze osservate per provare ad ottimizzare le stime dei
#   valori di aspettazione ottenuti

#   Allochiamo il post processor
    post_processor = POVMPostProcessor(result)

#   Settiamo il dual frame del post-porcessor
    post_processor.dual = \
        dual_from_empirical_frequencies(povm_post_processor=post_processor)

    return post_processor
    

def lancia_POVM(backend: BackendV2,         \
                circuit: QuantumCircuit,    \
                POVM: POVMImplementation)   \
                -> str:
#   Funzione che passato il circuito quantistico da implemntare e
#   l'implementazione delle POVM esegue la procedura desiderata

#   Creazione del PUB per il run della POVM
    pub = circuit

#   Inizializzazione del sampler per le POVM
    sampler = RuntimeSampler(mode=backend)

#   Attiviamo il dynamical decoupling con sequenza ...
    sampler.options.dynamical_decoupling.enable = True
    # sampler.options.dynamical_decoupling.sequence_type = "XX"
    sampler.options.dynamical_decoupling.sequence_type = "XpXm"
#     sampler.options.dynamical_decoupling.sequence_type = "XY4"

#   Allochiamo il POVM sampler
    POVM_sampler = POVMSampler(sampler=sampler)

#   Definiamo il numero di shots da eseguire
    shots = 20000

#   Eseguiamo circuito e POVM
    job = POVM_sampler.run(pubs=[pub], shots=shots, povm=POVM)

#   radice comune nomi file
    data = str(shots) + '_' + str(datetime.now())
    
#   Nome file contenente i metadata del job
    nome_file_metadata = './Jobs/POVM_metadata_shots_' + data + '.pkl'

#   Salviamo i metatdata del job --> il job_id è incluso nei metadata
    job.save_metadata(filename=nome_file_metadata)

#   Ritorniamo il la data univoca per recuperare job_id e relativi metadata
    return data

def recupera_POVM(service: QiskitRuntimeService, id: str) -> POVMPubResult:
#   Funzione che recupera dai server IBM i risultati del job relativo

#   Recuperiamo il job dai servere IBM
    nome_file = './Jobs/POVM_metadata_shots_' + id + '.pkl'
    job = POVMSamplerJob.recover_job(filename=nome_file, service=service)

#   Recuperiamo i risultati delle POVM
    result = job.result()

#   Ritorniamo i risultati scartando i metadate a noi non interessanti.
#   L'indice indica a quale pub ci riferiamo, nel nostro caso ne abbiamo
#   solo uno --> idx=0
    return result[0]