from qiskit.quantum_info import *
import numpy as np

#-----------------------------------------------------------------------------#
#   OK
#-----------------------------------------------------------------------------#

class HamOp:
    def __init__(self, n = 0):
        self.Hx_odd  = SparsePauliOp("I")
        self.Hy_odd  = SparsePauliOp("I")
        self.Hz_odd  = SparsePauliOp("I")

        self.Hx_even = SparsePauliOp("I")
        self.Hy_even = SparsePauliOp("I")
        self.Hz_even = SparsePauliOp("I")

        self.tot     = SparsePauliOp("I")
        
        self.n_spin  = n

    def genera_op_x(self, listaP_even, listaP_odd, pesi_even, pesi_odd):
        self.Hx_even = SparsePauliOp(listaP_even, pesi_even)
        self.Hx_odd  = SparsePauliOp(listaP_odd , pesi_odd )
        return

    def genera_op_y(self, listaP_even, listaP_odd, pesi_even, pesi_odd):
        self.Hy_even = SparsePauliOp(listaP_even, pesi_even)
        self.Hy_odd  = SparsePauliOp(listaP_odd , pesi_odd )
        return

    def genera_op_z(self, listaP_even, listaP_odd, pesi_even, pesi_odd):
        self.Hz_even = SparsePauliOp(listaP_even, pesi_even)
        self.Hz_odd  = SparsePauliOp(listaP_odd , pesi_odd )
        return
    
    def genera_op_tot(self, listaP, pesi):
        self.tot = SparsePauliOp(listaP, pesi)
        return

    def valuta_odd(self, J = (0,0,0), f = 1, param = 1, t = 0):
        return f(param, t) * (J[0] * self.Hx_odd.to_matrix()  \
                           +  J[1] * self.Hy_odd.to_matrix()) \
                           +  J[2] * self.Hz_odd.to_matrix()
        
    def valuta_even(self, J = (0,0,0), f = 1, param = 1, t = 0):
        return f(param, t) * (J[0] * self.Hx_even.to_matrix()  \
                           +  J[1] * self.Hy_even.to_matrix()) \
                           +  J[2] * self.Hz_even.to_matrix()
        
    def valuta_tot(self, J = (0,0,0), f = 1, param = 1, t = 0):
        return self.valuta_odd (J, f, param, t) \
             + self.valuta_even(J, f, param, t)

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

def crea_str_pauli(idx1, idx2, op, n):
#   Funzione che genera la stringa 'II...IIopII..IIopII..II'

#   Inizializzazione stringa vuota
    string = ''

#   Riempiamo la stringa nel modo opportuno
    for i in range(0,idx1):
        string += 'I'

    string += op

    for i in range(idx1+1,idx2):
        string += 'I'
    
    string += op

    for i in range(idx2+1,n):
        string += 'I'

    return string

#-----------------------------------------------------------------------------#

def crea_hamiltonianaOp(n, chiuso):
#   Funzione che crea le varie componenti dell'Hamiltoniana del sistema

#   Verifichiamo che la catena da costruire sia pari
    if n % 2 == 1:
        print('ERRORE:: Anello dispari')
        return 0, 0, 0, 0
    
#   Allochiamo l'oggetto contenente le varie componenti dell'Hamiltoniana
    ham = HamOp(n)

#   Inizializziamo le liste di Pauli dell'Hamiltoniana
    Px_even = []
    Py_even = []
    Pz_even = []
    Px_odd  = []
    Py_odd  = []
    Pz_odd  = []
    P_tot   = []

#   Cicliamo su tutte le stringhe generabili
    for i in range(0,n-1):

#       Creiamo gli operatori di Pauli opportuni.
#       L'operatore generato va diviso per 4 in quanto Sx = Pauli('X')/2 ect
        Px = crea_str_pauli(i, i+1, 'X', n)
        Py = crea_str_pauli(i, i+1, 'Y', n)
        Pz = crea_str_pauli(i, i+1, 'Z', n)

#       Sommiamo i vari operatori
        if i % 2 == 0:
            Px_even.append(Px)
            Py_even.append(Py)
            Pz_even.append(Pz)
        else:
            Px_odd.append(Px)
            Py_odd.append(Py)
            Pz_odd.append(Pz)
            
#   Verifichiamo se vogliamo un anello chiuso o aperto
    if chiuso == True:

#       Aggiungiamo il terine di interazione tra l'ultimo ed il primo sito
#       della catena [essendo n pari questo è sempre un termine dispari, 
#       visto che gli indici in python vanno da 0 a n-1]
        Px_odd.append((crea_str_pauli(0, n-1, 'X', n)))
        Py_odd.append((crea_str_pauli(0, n-1, 'Y', n)))
        Pz_odd.append((crea_str_pauli(0, n-1, 'Z', n)))

#   Inizializziamo gli opportuni pesi per ogli operatori da generare
#   Allochiamo i vettori dei pesi
    pesi_even = np.zeros(int(n/2)) + 1/4

    if chiuso == True:
        pesi_odd = np.zeros(int(n/2))     + 1/4
    else:
        pesi_odd = np.zeros(int(n/2) - 1) + 1/4

#   Generiamo gli opportuni operatori date le stringhe create
    ham.genera_op_x(Px_even, Px_odd, pesi_even, pesi_odd)
    ham.genera_op_y(Py_even, Py_odd, pesi_even, pesi_odd)
    ham.genera_op_z(Pz_even, Pz_odd, pesi_even, pesi_odd)

#   Concateniamo tutte le stringhe per poi generare la rappresentazione 
#   totale dell'Hamiltoniana
    P_tot = Px_even + Px_odd + Py_even + Py_odd + Pz_even + Pz_odd

#   Concateniamo tutti i pesi per avere un vettore delle dimensioni opportune
    pesi = np.concatenate([pesi_even,pesi_odd, pesi_even,pesi_odd, pesi_even,pesi_odd])

#   Generiamo la rappresentazione totale dell'Hamiltoniana
    ham.genera_op_tot(P_tot, pesi)

    return ham

#-----------------------------------------------------------------------------#

def group_potenza_hamiltoniana(Ham: HamOp, n: int) -> list[SparsePauliOp]:
#   Funzione che calcola la potenza n-esima dell'Hamiltoninana totale fornita.
#   Ad ogni moltiplicazione raggruppa gli elementi che commutano. In questo
#   modo la complessità del problema si riduce ma allo stesso tempo NON
#   otteniamo il raggruppamento ideale

#   Verificchiamo di dover calcolare una potenza effettiva
    if n == 0:
        return [SparsePauliOp(crea_str_pauli(1,2,'I',Ham.tot.num_qubits))]
    
#   Calcoliamo il primo grouping degli elementi che commutano di Ham
    group_H = Ham.tot.group_commuting()

#   Cicliamo sul numero di potenze da calcolare
    for _ in range(1,n):
        
#       Inizializziamo la lista di supporto
        group_tmp = []

        for j in range(0,len(group_H)):
#           Calcoliamo il prodotto tra l'elemento attuale e l'Hamiltoniana
            ham_tmp = group_H[j] @ Ham.tot

#           Aggiungiamo l'elemento totale all'Hamiltoniana
            group_tmp += ham_tmp.simplify().group_commuting()

#       Aggiorniamo il raggruppamento attuale
        group_H = group_tmp

#   Eseguiamo 

#   Ritorniamo il raggruppamento calcolato
    return group_H

#-----------------------------------------------------------------------------#

def costruisci_potenze_Ham(Ham: HamOp, pow: int) -> list[SparsePauliOp]:
#   Funzione che costruisce e potenze desiderate dell'Hamiltoniana

#   Allochiamo la lista che contrrà le potenze dell'Hamiltoniana
    Ham_pow = []

#   Aggiungiamo l'identità
    Ham_pow.append(SparsePauliOp(''.join(['I'] * Ham.n_spin)))

#   Controlliamo se dobbiamo fermarci o meno
    if pow < 1:
        return Ham_pow

#   Aggiungiamo l'Hamiltoniana alla lista
    Ham_pow.append(Ham)

#   Variabile d'appogio
    tmp = Ham.tot

#   Costruiamo tutte le altre potenze dell'Hamiltoniana
    for i in range(2, pow):
#       Costruiamo la potenza opportuna dell'Hamiltoniana
        tmp = tmp @ Ham.tot

#       Aggiungiamo la potenza alla lista degli operatori
        Ham_pow.append(tmp)
#   end i

    return Ham_pow

#-----------------------------------------------------------------------------#