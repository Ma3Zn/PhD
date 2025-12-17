module Circuit

include("Gates.jl")
using .Gates

include("Utils.jl")
using .Utils

export  apply_single_qubit_gate!, apply_two_qubit_gate!, make_adjacent!, 
        apply_heisenberg_2site!, apply_layer_heisenberg!, apply_trotter_II!, 
        apply_trotterized_annealing!

using ITensorMPS
using ITensors


"""
    apply_trotter_II!(mps, ctrl, targ, chi; cutoff)

Applica in-place (modificando `mps`) l'evoluzione dovuta all'intero processo di 
termalizzazione digitalizzato.
La mps è attesa già nell'opportuno stato iniziale del processo di annealing.
drive è la funzione di drive del processo adiabatico che trasforma il termine 
di Ising in un termine di Heisenberg.
drive_sp è una funzione che man mano va a spegne il campo presente sull'ultimo 
sito dell'anello per splittarne gli autovalori.
lam e lam_sp sono parametri che regolano la velocità dei due rispettivi drive.
`chi` è la dimensione massima del bond.
`cutoff` è il valore a cui troncare i valori singolari nella SVD.
"""
function apply_trotterized_annealing!(
    mps::MPS,               # MPS su cui applicare il circuito
    n_spin::Int,            # numero di qubit
    n_tstep::Int,           # numero di step Trotter
    J::Vector{<:Number},    # lista [Jx, Jy, Jz]
    drive::Function,        # funzione drive(param, t)
    lam::Float64,           # parametro lam per drive
    drive_sp::Function,     # funzione drive_sp(param_sp, n_spin, t)
    lam_sp::Float64,        # parametro lam per drive_sp
    chi::Int;               # bond dimension massima
    cutoff::Float64=1e-12
)
    # Se non sono richiesti trotter step ritorniamo senza apportare 
    # alcuna modifica
    if n_tstep == 0
        return
    end

    # 0) Parametri temporali dell'annealing
    #    tempo finale tf = 1/lam (se segui la formula del file python)
    tf = 1/lam
    dt = tf/n_tstep
    times = range(0, stop=tf-dt, length=n_tstep)

    # # 1) Stato iniziale: in Qiskit costruivi |1010...>. 
    # #    Qui presumo tu abbia già un MPS “|0000...>”
    # #    e vuoi replicare l'effetto di X sulle posizioni dispari:
    # for i in 1:2:n_spin-1
    #     apply_single_qubit_gate!(mps, i, X_gate(), chi)  # definisci X_gate() a 2x2

    # 2) Loop Trotter
    for step in 1:n_tstep
        t = times[step]
        t_eval = t + dt/2

        # Calcolo parametri:
        drive_val    = drive(lam, t_eval)                # es. drive(param, t)
        drive_sp_val = drive_sp(lam_sp, n_spin, t_eval)  # es. drive_sp(param_sp, n_spin, t)

        # Omelyan al II ordine:
        # crea_circuito_trotter_II(flag, dt, drive_val, drive_sp_val, J, n_spin)
        # In Qiskit facevi "flag = step%2" etc.
        flag = (step % 2)

        apply_trotter_II!(mps, flag, dt, drive_val, drive_sp_val, J, n_spin, chi; cutoff)
    end
end

#-----------------------------------------------------------------------------#

"""
    apply_trotter_II!(mps, ctrl, targ, chi; cutoff)

Applica in-place (modificando `mps`) l'evoluzione dovuta ad un intero step di 
trotterizzazione tra siti pari e siti dispari dell'hamiltoniana.
In particolare il valore fornito in flag individua se eseguire una trotterizzazione 
partendo da siti pari o siti dispari.
Questo sarà passato 0 per i siti "dispari" quindi primo ctlr in 1, mentre sarà 
passato 1 per i siti "pari" quindi primo ctrl in 2 --- maledetto Julia e tutti
altri linguaggi che iniziano ad incizzare da 1 e non da 0 .........
`chi` è la dimensione massima del bond.
`cutoff` è il valore a cui troncare i valori singolari nella SVD.
"""
function apply_trotter_II!(
    mps::MPS,
    flag::Int,
    dt::Float64, 
    drive_val::Float64, 
    drive_sp_val::Float64,
    J::Vector{<:Number}, 
    n_spin::Int,
    chi::Int;
    cutoff::Float64=1e-12
)
    # step 0: singolo corpo sull'ultimo qubit (Rz)
    if drive_sp_val > 0
        # Rz(drive_sp_val * 0.6 * dt/2)
        # ATTENZIONE::  lo 0.6 viene da come definiamo l'hamiltoniana a 
        #               singolo corpo in MATLAB
        gate = rz_gate(drive_sp_val * 0.6 * dt/2)
        apply_single_qubit_gate!(mps, n_spin, gate)
    end

    # step 1
    apply_layer_heisenberg!(mps, flag,  0.1931833275037836*dt, drive_val, J, n_spin, chi; cutoff)

    # step 2
    apply_layer_heisenberg!(mps, flag+1, 0.5*dt, drive_val, J, n_spin, chi; cutoff)

    # step 3
    apply_layer_heisenberg!(mps, flag,  0.613633344992433*dt, drive_val, J, n_spin, chi; cutoff)

    # step 4
    apply_layer_heisenberg!(mps, flag+1, 0.5*dt, drive_val, J, n_spin, chi; cutoff)

    # step 5
    apply_layer_heisenberg!(mps, flag,  0.1931833275037836*dt, drive_val, J, n_spin, chi; cutoff)

    # step 6: di nuovo singolo corpo
    if drive_sp_val > 0
        gate = rz_gate(drive_sp_val * 0.6 * dt/2)
        apply_single_qubit_gate!(mps, n_spin, gate)
    end
end

#-----------------------------------------------------------------------------#

"""
    apply_layer_heisenberg!(mps, ctrl, targ, chi; cutoff)

Applica in-place (modificando `mps`) l'evoluzione dovuta ad un intero layer di
termini di Heisenberg su i siti pari/dispari di tutto l'anello.
In particolare la parità/disparità dei sti viene individuata dal valore di flag.
`chi` è la dimensione massima del bond.
`cutoff` è il valore a cui troncare i valori singolari nella SVD.
"""
function apply_layer_heisenberg!(
    mps::MPS,
    flag::Int,
    dt_layer::Float64,
    drive_val::Float64,
    J::Vector{<:Number},    # [Jx, Jy, Jz]
    n_spin::Int,
    chi::Int;
    cutoff::Float64=1e-12
)
    # Rinormalizziamo il valore del flag.
    # Questo perchè allo step precedente se partiamo da un indice dispari (0)
    # verrà passato flag = 0/1. Mentre se partiamo da un indice pari (1) verrà
    # passato un flag = 1/2.
    # Maledetto chi ha pensato di far partire l'indicizzazione da 1 ..........
    flag = flag % 2

    # In Qiskit: J_tmp[0] = J[0]*drive_val, 
    #            J_tmp[1] = J[1]*drive_val,
    #            J_tmp[2] = J[2]
    # Prepara i coefficienti "temporanei" J_tmp
    J_tmp = [J[1]*drive_val, J[2]*drive_val, J[3]] 

    # Cicliamo su tutte le coppie di siti pari/dispari, in base al valore fornito di flag
    # Il + 1 è necessario in quanto flag = 0/1
    i = flag + 1
    while i < n_spin

        # Se ci troviamo sul bordo e la catena è "chiusa"
        if i == n_spin
            j = 1
        else
            j = i + 1
        end

        apply_heisenberg_2site!(mps, i, j, dt_layer, J_tmp, chi; cutoff)

        i += 2
    end
end

#-----------------------------------------------------------------------------#

"""
    apply_heisenberg_2site!(mps, ctrl, targ, chi; cutoff)

Applica in-place (modificando `mps`) l'evoluzione dovuta al termine di 
Heisenberg attuale tra il sito crtl e targ.
In particolare avremo sempre che targ = ctrl + 1 tranne quando ctrl = N 
e targ = 1, i.e., quando siamo al bordo periodico dell'anello.
`chi` è la dimensione massima del bond.
`cutoff` è il valore a cui troncare i valori singolari nella SVD.
"""
function apply_heisenberg_2site!(
    mps::MPS,
    ctrl::Int,
    targ::Int,
    dt::Float64,
    J::Vector{<:Number}, # (Jx, Jy, Jz)
    chi::Int;
    cutoff::Float64=1e-12
)
    # 1) Verifichiamo se i due siti sono adiacenti o se dobbiamo implementare la
    #    condizione periodica al bordo
    if targ - ctrl == 1

        # Calcoliamo l'opportuna rappresentazione matriciale dei CNOT da usare
        #
        #   ctrl --------*--------
        #                |
        #                |
        #   targ --------X--------
        CNOT_rep = cnot_gate_clt()
        
        # Inizializziamo il flag che ci dice se al termine del layr dobbiamo o
        # meno riportare in prima posizione il ctarget attuale
        scambio = false
    else
        # Riordiniamo in modo opportuno control e target
        tmp  = targ
        targ = ctrl
        ctrl = tmp

        # Rendiamo adiacenti i due siti
        make_adjacent!(mps, ctrl, target, chi; cutoff)
        
        # Aggiorniamo in modo opportuno l'indice del ctrl
        ctrl = target - 1

        # Calcoliamo l'opportuna rappresentazione matriciale dei CNOT da usare
        #
        #   targ --------X--------
        #                |
        #                |
        #   ctrl --------*--------
        CNOT_rep = cnot_gate_cgt()
        
        # Inizializziamo il flag che ci dice se al termine del layr dobbiamo o
        # meno riportare in prima posizione il ctarget attuale
        scambio = true
    end


    # 2) Applica gate "CNOT(ctrl, targ)"
    apply_two_qubit_gate!(mps, ctrl, CNOT_rep, chi; cutoff)

    # 3) rx((Jx*dt - pi)/2, ctrl)
    apply_single_qubit_gate!(mps, ctrl, rx_gate( (J[1]*dt - π) / 2) )
    # 4) h(ctrl)
    apply_single_qubit_gate!(mps, ctrl, h_gate())
    # 5) rz(Jz*dt/2, targ)
    apply_single_qubit_gate!(mps, targ, rz_gate( J[3]*dt / 2 ))

    # 6) cx(ctrl, targ)
    apply_two_qubit_gate!(mps, ctrl, CNOT_rep, chi; cutoff)

    # 7) h(ctrl)
    apply_single_qubit_gate!(mps, ctrl, h_gate())
    # 8) rz(-Jy*dt/2, targ)
    apply_single_qubit_gate!(mps, targ, rz_gate( -J[2]*dt / 2 ))

    # 9) cx(ctrl, targ)
    apply_two_qubit_gate!(mps, ctrl, CNOT_rep, chi)

    # 10) rx(pi/2, ctrl)
    apply_single_qubit_gate!(mps, ctrl, rx_gate(π/2))
    # 11) rx(-pi/2, targ)
    apply_single_qubit_gate!(mps, targ, rx_gate(-π/2))

    # 12) Controlliamo se dobbiamo riportare il target al primo sito
    # dell'anello
    if scambio
        # Riportiamo il target attuale nel primo sito dell'anello
        # che attualmente è nella variabile ctrl perchè invochiamo il
        # CNOT nella variante cgt (ovvero quella in cui il control 
        # effettivo è in posizione ctrl+1)
        make_adjacent!(mps, ctrl, 1, chi; cutoff)
    end
end

#-----------------------------------------------------------------------------#

"""
    make_adjacent!(mps, i, j, chi; cutoff)

Applica in-place (modificando `mps`) una catena di swap per portare il sito i 
in posizione j-1.
`chi` è la dimensione massima del bond.
`cutoff` è il valore a cui troncare i valori singolari nella SVD.
"""
function make_adjacent!(
    mps::MPS,
    i::Int,
    j::Int,
    chi::Int;
    cutoff::Float64=1e-12
)

    # Verifichiamo lordinamento degli SWAP da eseguire, 
    # ovvero se dobbiamo avvicinare da sinistra o da destra i a j
    if i < j
        # Avviciniamo da destra i a j.
        # In particolare i = 1 e j = N --> i = N-1 e j = N
        for k in i:j-2
            apply_two_qubit_gate!(mps, k, swap_gate(), chi; cutoff)
        end
    else
        # Avviciniamo da sinistra i a j. In particolare per i nostri usi 
        # (ovvero riportare in posizione 1 il primo qubit) vogliamo eseguire
        # uno swap effettivo anche tra i e j in questo caso.
        # In termini più semplici j = 1 e i = N-1 --> i = 1 e j = 2
        for k in reverse(j:i-1)
            apply_two_qubit_gate!(mps, k, swap_gate(), chi; cutoff)
        end
    end

    # ATTENZIONE:: Essendo una funzione in-place credo che questo return sia 
    #              inutile
    return mps
end

#-----------------------------------------------------------------------------#

"""
    apply_single_qubit_gate!(mps, i, gate_matrix, chi)

Applica in-place (modificando `mps`) un gate 2×2 al sito `i`.
`gate_matrix` è una Matrix{<:Real/Complex} 2×2.
"""
function apply_single_qubit_gate!(
    mps::MPS,
    i::Int,
    gate_matrix::AbstractMatrix{<:Number}
)
    # 1) Creare l'ITensor corrispondente al gate (indice in, indice out).   
    s_in  = siteind(mps, i) # indice fisico "in"
    s_out = prime(s_in)     # indice fisico "out" (versione primata)
    GateT = ITensor(gate_matrix, s_out, s_in)

    # 2) Contraiamo i tensori e deprimiamo gli indici
    mps[i] = noprime(mps[i] * GateT)

    # ATTENZIONE:: Essendo una funzione in-place credo che questo return sia 
    #              inutile
    return mps
end

#-----------------------------------------------------------------------------#

"""
apply_two_qubit_gate!(
    mps::MPS,
    i::Int,
    gate_matrix::AbstractMatrix{<:Number},
    chi::Int;
    cutoff=1e-12
)

Applica un gate 4×4 ai siti i e i+1 di un MPS (qubit).
Dopo la contrazione, suddivide manualmente gli indici,
in modo che A_i e A_{i+1} abbiano ciascuno UN indice fisico (dim=2).
Funziona anche nel caso di 2 siti totali con bond=1 (nessun link).
"""
function apply_two_qubit_gate!(
    mps::MPS,
    i::Int,
    gate_matrix::AbstractMatrix{<:Number},
    chi::Int;
    cutoff::Float64=1e-12
)
    # To apply this gate in a controlled manner, first 'gauge' the 
    # MPS psi such that either site i or i+1 is the orthogonality center.
    orthogonalize!(mps, i)

    # Recuperiamo un ITensor dalla rappresentazione matriciale del nostro gate
    # 2) Gate 4×4 come ITensor:
    s_i   = siteind(mps, i)      # dimension 2
    s_ip1 = siteind(mps, i+1)    # dimension 2
    s_i_out   = prime(s_i)
    s_ip1_out = prime(s_ip1)
    GateT = ITensor(gate_matrix, s_ip1, s_i, s_ip1_out, s_i_out)

    # The other MPS tensors are now either left-orthogonal or right-orthogonal 
    # and can be left out of further steps without producing incorrect results.
    #
    # Next, contract the gate tensor G with the MPS tensors
    contraction = (mps[i] * mps[i+1]) * GateT
    contraction = noprime(contraction)
    
    # Finally, use the singular value decomposition (SVD) to factorize the 
    # resulting tensor, multiplying the singular values into either U or V.
    # Assign these two tensors back into the MPS to update it.
    inds = uniqueinds(mps[i], mps[i+1])
    U,S,V = svd(contraction, inds; cutoff=cutoff, maxdim=chi)
    mps[i]   = U
    mps[i+1] = S*V

    # ATTENZIONE:: Essendo una funzione in-place credo che questo return sia 
    #              inutile
    return mps
end

end # module

