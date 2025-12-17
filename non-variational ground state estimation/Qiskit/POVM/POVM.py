from qiskit import *
from qiskit.primitives import *
from qiskit.quantum_info import *

from povm_toolbox.library import *
from povm_toolbox.sampler import *
from povm_toolbox.post_processor import *

from qiskit.primitives import StatevectorSampler

def alloca_POVM_post_processor(result: POVMPubResult) -> POVMPostProcessor:
#   Funzione che alloca il POVMPostProcessor necessario per eseguire la stima
#   dei valori di aspettazione delle osservabili una volta che abbiamo eseguito
#   le POVM. In particolare setiamo il dual frame del post-processor sulla
#   strategia delle frequenze osservate per provare ad ottimizzare le stime dei
#   valori di aspettazione ottenuti

#   Allochiamo il post processor
    post_processor = POVMPostProcessor(result)

#   Settiamo il dual frame del post-porcessor
    # post_processor.dual = \
    #     dual_from_empirical_frequencies(povm_post_processor=post_processor)

    return post_processor
    

def esegui_POVM(circuit: QuantumCircuit, POVM: POVMImplementation) \
                -> POVMPubResult:
#   Funzione che passato il circuito quantistico da implemntare e
#   l'implementazione delle POVM esegue la procedura desiderata

#   Creazione del PUB per il run della POVM
    pub = circuit

#   Inizializzazione del sampler per le POVM
    sampler = StatevectorSampler()
    POVM_sampler = POVMSampler(sampler=sampler)

#   Definiamo il numero di shots da eseguire
    shots = 10000

#   Eseguiamo circuito e POVM
    job = POVM_sampler.run(pubs=[pub], shots=shots, povm=POVM)

#   Recuperiamo i risultati delle POVM
    result = job.result()

#   Ritorniamo. i risultati scartando i metadate a noi non interessanti.
#   L'indice indica a quale pub ci riferiamo, nel nostro caso ne abbiamo
#   solo uno --> idx=0
    return result[0]