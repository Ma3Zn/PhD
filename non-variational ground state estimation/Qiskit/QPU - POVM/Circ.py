import numpy as np

from qiskit import *

from math import pi

#-----------------------------------------------------------------------------#
#   OK
#-----------------------------------------------------------------------------#

def crea_circuito_inizializzazione(n):
#   Funzione che genera il circuito per preparare lo stato |1010.10> -- dove 
#   si è usata la convenzione di Qiskit di labellare i qubit da destra a sinistra, 
#   i.e. il qubit 0 (il primo dei pari) è quello che varia più velocemente

#   Allochiamo il nostro circuito
    qc = QuantumCircuit(n,n)

#   Dobbiamo applicare un X a tutti i qubi in posizione dispari
    id_q = range(1,n,2)
    
#   Eseguiamo lo scheduling dei gate
    for i in id_q:
        qc.x(i)

#   Ritorniamo il circuito generato
    return qc

#-----------------------------------------------------------------------------#

def crea_termine_Heisenberg(ctrl, target, dt, J, n_spin):
#   Funzione che genera il circuito per l'implemntazione di un'Hamiltoniana
#   con interazione di Heisenberg

#   Allochiamo il nostro circuito
    qc = QuantumCircuit(n_spin, n_spin)

#   Verifichiamo che il target faccia parte del circuito e nel caso
#   implementiamo la chiusura periodica dell'anello
    if target > (n_spin - 1):
        target = 0

#   Implementiamo il circuito peril termine di Heisenberg usando la 
#   decomposizione minimale in 3-CNOT che puoi trovare sulla tesi di dottorato 
#   di Tacchino oppure in "G. Vidal and C. M. Dawson, Universal quantum circuit
#   for two-qubit transformations with three controlled-NOT gates,
#   Phys. Rev. A 69, 010301 (2004)."
    
    qc.cx(ctrl, target)

    qc.rx((J[0]*dt - pi) / 2, ctrl)
    qc.h(ctrl)
    qc.rz(J[2]*dt / 2, target)

    qc.cx(ctrl, target)

    qc.h(ctrl)
    qc.rz(-J[1]*dt / 2, target)

    qc.cx(ctrl, target)

    qc.rx(pi/2, ctrl)
    qc.rx(-pi/2, target)

#   Ritorniamo il circuito generato
    return qc

#-----------------------------------------------------------------------------#

def crea_layer_Heisenberg(flag, dt, drive_val, J, n_spin, t_coeff):
#   Funzione che genera un layer per interazioni Heisenberg su tutti
#   i qubit individuati dal flag passato

#   Allochiamo il nostro circuito
    qc = QuantumCircuit(n_spin, n_spin)

#   Rinormalizziamo il flag
    flag = flag % 2

#   Calcoliamo il valore opportuno di Jx [J[0]] e Jy [J[1]]
    J_tmp = [0,0,0]
    J_tmp[0] = J[0] * drive_val
    J_tmp[1] = J[1] * drive_val
    J_tmp[2] = J[2]

#   Calcoliamo il valore opportuno di dt
    dt_tmp = dt * t_coeff

#   Generiamo tutte le coppie di qudit
    for id in range(flag, n_spin, 2):

#       generiamo il circuito relativo al termine attuale
        qc.compose(crea_termine_Heisenberg(id, id + 1, dt_tmp, J_tmp, n_spin), inplace=True)

#   Ritorniamo il circuito generato
    return qc

#-----------------------------------------------------------------------------#

def crea_circuito_trotter_II(flag, dt, drive_val, drive_sp_val, J, n_spin):
#   Funzione che genera il circuito per un passo dell'annealing digitalizzato
#   trotterizzando a sua volta le interazioni pari e quelle dispari tramite la
#   formula al second'ordine di Omelyan

#   Allochiamo il nostro circuito
    qc = QuantumCircuit(n_spin, n_spin)

#   Verifichiamo che il termine a singolo corpo sia effettivamente presente
    if drive_sp_val > 0:
#       Trotterizziamo il passo a singoo corpo sull'ultimo qubit
        qc.rz(drive_sp_val * 0.6 * dt/2, n_spin-1)

#   Generiamo il primo step della formula di omelyan su tutti i qubit
#   indicati dal flag passato
    qc.compose(crea_layer_Heisenberg(flag, dt, drive_val, J, n_spin, 0.1931833275037836), inplace=True)

#   Generiamo il secondo step della formula di omelyan su tutti i qubit
#   indicati dal flag passato
    qc.compose(crea_layer_Heisenberg(flag + 1, dt, drive_val, J, n_spin, 0.5), inplace=True)

#   Generiamo il terzo step della formula di omelyan su tutti i qubit
#   indicati dal flag passato
    qc.compose(crea_layer_Heisenberg(flag, dt, drive_val, J, n_spin, 0.613633344992433), inplace=True)

#   Generiamo il quarto step della formula di omelyan su tutti i qubit
#   indicati dal flag passato
    qc.compose(crea_layer_Heisenberg(flag + 1, dt, drive_val, J, n_spin, 0.5), inplace=True)

#   Generiamo il quinto step della formula di omelyan su tutti i qubit
#   indicati dal flag passato
    qc.compose(crea_layer_Heisenberg(flag, dt, drive_val, J, n_spin, 0.1931833275037836), inplace=True)

#   Verifichiamo che il termine a singolo corpo sia effettivamente presente
    if drive_sp_val > 0:
#       Trotterizziamo il passo a singoo corpo sull'ultimo qubit
        qc.rz(drive_sp_val * 0.6 * dt/2, n_spin-1)

#   Ritorniamo il circuito generato
    return qc
    
#-----------------------------------------------------------------------------#

def crea_circuito_annealing(J, n_spin, n_tstep, f, param, f_sp, param_sp):
#   Funzione che crea il circuito per la digitalizzazione
#   [in n_tstep trotter step] del processo di annealing

#   Allochiamo il nostro circuito
    qc = QuantumCircuit(n_spin,n_spin)

    if n_tstep == 0:
        return qc

#   Calcoliamo il tempo finale dell'evoluzione influenzato dal parametro
#   di velocità dell'annealing
    tf = 1/param

#   Calcoliamo il dt di ogni trotter step
    dt = tf / n_tstep

#   Calcoliamo i vari istanti "iniziali" di ogni trotter step
    time = np.linspace(0, tf - dt, n_tstep)

#   Cicliamo su tutti i trotter step da eseguire
    for it in range(0,time.size):

#       Recuperiamo il tempo attuale
        t = time[it]

#       Calcoliamo il centro di valutazione temporale di ogni step
        t_eval = t + dt/2

#       Calcoliamo il valore del drive dell'annealing per il trotter step
        drive_val = f(param, t_eval)

#       Calcoliamo il valore del drive dell'annealing per il termine a
#       singolo corpo applicato al solo ultimo qubit della catena
        drive_sp_val = f_sp(param_sp, n_spin, t_eval)

#       Per il trotter step attuale creaimo il circuito per la 
#       trotterizzazione di Omelyan al II ordine per H_odd ed H_even.
#       Come flag passiamo it%2 + 1 per essere coerenti con MATLAB per
#       poter verificare i risultati
#       TODO:: Una volta verificati i risultati con MATLAB leva questo +1
#       per aumentare la coerenza dei vari script python
        qc.compose(crea_circuito_trotter_II(it % 2, dt, drive_val, drive_sp_val, J, n_spin), inplace=True)

#   Ritorniamo il circuito generato
    return qc

#-----------------------------------------------------------------------------#