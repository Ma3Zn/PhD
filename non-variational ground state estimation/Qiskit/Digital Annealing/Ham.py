from qiskit.quantum_info import Pauli
import numpy as np
#-----------------------------------------------------------------------------#
#   OK
#-----------------------------------------------------------------------------#

class Ham:
    def __init__(self, dim = 0):
        self.Hx_odd  = np.matrix(np.zeros((dim,dim))) * 1j
        self.Hy_odd  = np.matrix(np.zeros((dim,dim))) * 1j
        self.Hz_odd  = np.matrix(np.zeros((dim,dim))) * 1j

        self.Hx_even = np.matrix(np.zeros((dim,dim))) * 1j
        self.Hy_even = np.matrix(np.zeros((dim,dim))) * 1j
        self.Hz_even = np.matrix(np.zeros((dim,dim))) * 1j

    def valuta_odd(self, J = (0,0,0), f = 1, param = 1, t = 0):
        return f(param, t) * (J[0] * self.Hx_odd + J[1] * self.Hy_odd) \
               +              J[2] * self.Hz_odd
        
    def valuta_even(self, J = (0,0,0), f = 1, param = 1, t = 0):
        return f(param, t) * (J[0] * self.Hx_even + J[1] * self.Hy_even) \
               +              J[2] * self.Hz_even
        
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

def crea_hamiltoniana(n, chiuso):
#   Funzione che crea le varie componenti dell'Hamiltoniana del sistema

#   Verifichiamo che la catena da costruire sia pari
    if n % 2 == 1:
        print('ERRORE:: Anello dispari')
        return 0, 0, 0, 0
    
#   Calcoliamo la dimensione del sistema finale
    dim = 2**n
    
#   Allochiamo l'oggetto contenente le varie componenti dell'Hamiltoniana
    ham = Ham(dim)

#   Cicliamo su tutte le stringhe generabili
    for i in range(0,n-1):

#       Creiamo gli operatori di Pauli opportuni.
#       L'operatore generato va diviso per 4 in quanto Sx = Pauli('X')/2 ect
        Px = np.matrix(Pauli(crea_str_pauli(i, i+1, 'X', n)).to_matrix() / 4)
        Py = np.matrix(Pauli(crea_str_pauli(i, i+1, 'Y', n)).to_matrix() / 4)
        Pz = np.matrix(Pauli(crea_str_pauli(i, i+1, 'Z', n)).to_matrix() / 4)

#       Sommiamo i vari operatori
        if i % 2 == 0:
            ham.Hx_even += Px
            ham.Hy_even += Py
            ham.Hz_even += Pz
        else:
            ham.Hx_odd  += Px
            ham.Hy_odd  += Py
            ham.Hz_odd  += Pz
            
#   Verifichiamo se vogliamo un anello chiuso o aperto
    if chiuso == True:

#       Aggiungiamo il terine di interazione tra l'ultimo ed il primo sito
#       della catena [essendo n pari questo è sempre un termine dispari, 
#       visto che gli indici in python vanno da 0 a n-1]
        ham.Hx_odd += \
            np.matrix(Pauli(crea_str_pauli(0, n-1, 'X', n)).to_matrix() / 4)
        ham.Hy_odd += \
            np.matrix(Pauli(crea_str_pauli(0, n-1, 'Y', n)).to_matrix() / 4)
        ham.Hz_odd += \
            np.matrix(Pauli(crea_str_pauli(0, n-1, 'Z', n)).to_matrix() / 4)

    return ham

#-----------------------------------------------------------------------------#