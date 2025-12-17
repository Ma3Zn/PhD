#!/usr/bin/env julia
using Test
using ITensors
using ITensorMPS

# --------------------------------------------------------------------------
# Includi i file con le definizioni di gate e funzioni di applicazione:
# --------------------------------------------------------------------------
include("../src/Circuit.jl")   # se definisce apply_single_qubit_gate!, apply_two_qubit_gate!, ecc.
using .Circuit

include("../src/Gates.jl")     # definisce x_gate, y_gate, z_gate, h_gate, rx_gate, rz_gate, swap_gate, ...
using .Gates

# Se hai già definito i gate two-qubit (CNOT, iSWAP, CZ), bene; altrimenti puoi inserirli qui.
# Se hai definito "States.jl" con varie funzioni, puoi includerle. Altrimenti definisci tutto qui.

###############################################################################
# Funzioni di utilità per questi test
###############################################################################

"""
    create_product_mps(bitstring) -> MPS

Crea un MPS (con siteinds("Qubit", length(bitstring))) tale che ogni qubit
è |0> o |1> a seconda del carattere in `bitstring`.
Esempio: bitstring="010" => |0,1,0>.
"""
function create_product_mps(bitstring::String)
    N = length(bitstring)
    sites = siteinds("S=1/2", N)
    # productMPS ci consente di specificare uno stato di prodotto
    # con un array di stringhe corrispondenti ai singoli siti.
    states = [c == '0' ? "Up" : "Dn" for c in bitstring]
    return MPS(sites, states)
end

"""
    mps_to_vector(mps::MPS) -> Vector{ComplexF64}

Estrae i coefficienti di `mps` (che ha N qubit) in un vettore
di dimensione 2^N, dove l'indice k (0-based) è associato
allo stato |b1 b2 ... bN> con b1 come bit più significativo
e bN meno significativo.

Ordine degli indici fisici in ITensor:
   - site(mps,1) = qubit 1
        # Per chiarezza: definisco le 8 stringhe possibili
        const ALL_BITSTRINGS_3 = ("000","001","010","011","100","101","110","111")
   - site(mps,2) = qubit 2
   ...
   - site(mps,N) = qubit N
"""
function mps_to_vector(mps::MPS)
    N = length(mps)
    # Contraiamo tutti i tensori dell'MPS in uno solo
    bigpsi = mps[1]
    for i in 2:N
        bigpsi *= mps[i]
    end

    # Preallochiamo il vettore di dimensione 2^N
    dim = 2^N
    psi_vec = Vector{ComplexF64}(undef, dim)

    # Ricaviamo gli indici fisici [site(mps,1), ..., site(mps,N)] nell'ordine
    s = [siteind(mps, i) for i in 1:N]  # oppure siteinds(mps)[i], dipende dalla versione ITensors

    # Funzione helper per estrarre il bit i-esimo (0=LSB, N-1=MSB)
    getbit(k, pos) = (k >> pos) & 1

    for k in 0:(dim-1)
        # vogliamo b1 = bit più significativo => bit(N-1)
        # b2 => bit(N-2), ...
        # bN => bit(0)
        # Poi passiamo (s1 => b1+1, s2 => b2+1, ..., sN => bN+1)
        # ad bigpsi[...] per ottenere l'ampiezza.

        # costruiamo la tupla di indici (b1+1, b2+1, ..., bN+1)
        cinds = ntuple(i -> getbit(k, N - i) + 1, N)
        #  i=1 => bit(N-1)
        #  i=2 => bit(N-2)
        #  ...
        #  i=N => bit(0)

        psi_vec[k+1] = bigpsi[cinds...]
    end

    return psi_vec
end

"""
    check_mps_state(mps, expected::Vector{ComplexF64}; atol=1e-10)

Confronta la wavefunction codificata in `mps` con quella fornita da `expected`.
Se la norma della differenza è < atol, il test passa.
"""
function check_mps_state(mps::MPS, expected::Vector{<:Number}; atol=1e-10)
    actual = mps_to_vector(mps)
    @test norm(actual - expected) < atol
end

"""
    bitstring_to_index3(bs::String) -> Int

Converte una stringa bs di lunghezza 3 (es. "010") in un indice [1..8],
interpretando bs[1] come bit più significativo, e bs[3] come bit meno significativo.
Quindi "000" => 1, "001" => 2, "010" => 3, "011" => 4,
     "100" => 5, "101" => 6, "110" => 7, "111" => 8.
"""
function bitstring_to_index3(bs::String)
    @assert length(bs) == 3 "La stringa deve essere di lunghezza 3"
    # estrai i bit come 0 o 1
    a = (bs[1] == '1') ? 1 : 0
    b = (bs[2] == '1') ? 1 : 0
    c = (bs[3] == '1') ? 1 : 0
    # calcola l'indice in [1..8]
    # a e' il bit piu' significativo
    # l'indice si ottiene da (4*a + 2*b + c + 1)
    return 4*a + 2*b + c + 1
end

"""
    build_expected_for_rxrz_3qubit(gkind, qubit, θ, initbits) -> Vector{ComplexF64}

Data una stringa di 3 bit (initbits), costruisce il vettore 8D corrispondente,
poi applica un gate single-qubit:
 - if gkind=="Rx", usa rx_gate(θ)
 - if gkind=="Rz", usa rz_gate(θ)
sul qubit in posizione `qubit` ∈ {1,2,3}, restituendo il vettore finale.
"""
function build_expected_for_rxrz_3qubit(gkind::String, qubit::Int, θ::Real, initbits::String)
    # 1) Costruisce lo stato iniziale in forma di vettore 8D
    vec_init = build_3qubit_vector(initbits)

    # 2) Seleziona la matrice 2×2 del gate
    local mat
    if gkind == "Rx"
        mat = rx_gate(θ)
    elseif gkind == "Rz"
        mat = rz_gate(θ)
    else
        error("build_expected_for_rxrz_3qubit: gkind='$gkind' non supportato (usa 'Rx' o 'Rz')")
    end

    # 3) Applica la single-qubit op globale (es: mat ⊗ I ⊗ I, se qubit=1, ecc.)
    vec_final = apply_single_qubit_op_global(mat, qubit, vec_init)

    return vec_final
end

"""
    build_3qubit_vector(bitstring::String) -> Vector{ComplexF64}

Ritorna un vettore di dimensione 8 (2^3) con 1 nella posizione corrispondente a `bitstring`
e 0 altrove. Interpreta bitstring[1] come bit più significativo,
bitstring[3] come bit meno significativo.
Esempio: "000" => [1,0,0,0,0,0,0,0]^T, "111" => [0,0,0,0,0,0,0,1]^T, ecc.
"""
function build_3qubit_vector(bitstring::String)::Vector{ComplexF64}
    @assert length(bitstring) == 3 "bitstring deve essere di lunghezza 3"

    # Converte la bitstring (es. "010") in un intero idx in [1..8]
    idx = bitstring_to_index3(bitstring)

    # Alloca vettore 8D
    vec = zeros(ComplexF64, 8)

    # Mettiamo 1 in posizione idx
    vec[idx] = 1 + 0im

    return vec
end

"""
Genera l'identita' 2x2
"""
function I2()
    return [1 0;
            0 1]
end

"""
    apply_single_qubit_op_global(mat2x2, qubitPos, vec8) -> Vector{ComplexF64}

Applica la matrice 2×2 `mat2x2` al qubit in posizione `qubitPos` ∈ {1,2,3}
di uno stato globale a 3 qubit, rappresentato da `vec8` di dimensione 8.

Le posizioni qubitPos=1,2,3 corrispondono a:
  - qubit1: bit più significativo
  - qubit2: bit intermedio
  - qubit3: bit meno significativo

Restituisce il nuovo stato globale come Vector{ComplexF64} di dimensione 8.
"""
function apply_single_qubit_op_global(mat2x2::AbstractMatrix, qubitPos::Int, vec8::Vector{ComplexF64})
    # Costruiamo l'operatore 8×8
    # se qubitPos=1 => mat2x2 ⊗ I2 ⊗ I2
    # se qubitPos=2 => I2 ⊗ mat2x2 ⊗ I2
    # se qubitPos=3 => I2 ⊗ I2 ⊗ mat2x2
    Id = I2()  # identità 2×2
    global_op = nothing

    if qubitPos == 1
        # qubit1 => mat2x2, qubit2 => Id, qubit3 => Id
        global_op = kron3(mat2x2, Id, Id)
    elseif qubitPos == 2
        global_op = kron3(Id, mat2x2, Id)
    elseif qubitPos == 3
        global_op = kron3(Id, Id, mat2x2)
    else
        error("apply_single_qubit_op_global: qubitPos deve essere in {1,2,3}")
    end

    @assert size(global_op,1)==8 && length(vec8)==8 "dimensione mismatch"

    # Ora moltiplichiamo
    newvec = global_op * vec8
    return newvec
end

"""
    kron3(A, B, C) -> Matrix

Calcola il prodotto di Kronecker A x B x C di tre matrici.
"""
function kron3(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix)
    return kron(A, kron(B, C))
end

# ----------------------------------------------------------------------
    # 2) Definiamo le utilità per costruire un operatore 8×8
    #    che applichi "mat2q (4×4)" su (1,2) o (2,3)
    # ----------------------------------------------------------------------

    """
        build_2qubit_global_op(mat2q::AbstractMatrix, pair::Tuple{Int,Int}) -> 8×8

    Costruisce l'operatore 8×8 che agisce come mat2q (4×4) su qubit (pair),
    e come identita` sugli altri.
    - pair=(1,2) => op = (mat2q)⊗I2
    - pair=(2,3) => op = I2⊗(mat2q)
    Qui assumiamo 3 qubit => dimensione 2^3=8.
    """
    function build_2qubit_global_op(mat2q, pair::Tuple{Int,Int})
        # qubit1 => "left", qubit2 => "middle", qubit3 => "right"
        # dimension 2 each => total 8.
        if pair==(1,2)
            # mat2q su (q1,q2), I2 su q3 => op = kron3(mat2q, I2())
            # cioe` mat2q e` 4×4 => eq. "2 qubit", e q3 => dimension 2
            return kron(mat2q, I2())
        elseif pair==(2,3)
            # op = kron3(I2(), mat2q)
            return kron(I2(), mat2q)
        else
            @error "Pair $(pair) non gestito!"
        end
    end

    """
        apply_2qubit_op_global(mat2q, pair, initvec8) -> finalvec8

    Applica l'operatore 8×8 corrispondente a mat2q su (pair),
    partendo da initvec8, restituendo finalvec8.
    """
    function apply_2qubit_op_global(mat2q, pair, initvec8)
        op = build_2qubit_global_op(mat2q, pair)
        return op*initvec8
    end

    function build_3qubit_vector(bs::String)
        # produce un Vector{ComplexF64} di length=8 con 1 in posizione corrispondente, 0 altrove.
        vec = zeros(ComplexF64,8)
        idx = bitstring_to_index3(bs) # come in precedenza => 1..8
        vec[idx] = 1
        return vec
    end

# --------------------------------------------------------------------------
# 2) Scelta di bitstring da testare
# --------------------------------------------------------------------------
# Se vuoi testare TUTTI gli 8 stati, li elenchi:
const ALL_BITSTRINGS_3 = ("000","001","010","011","100","101","110","111")

# --------------------------------------------------------------------------
# 3) Test Single-Qubit gates (X, Y, Z, H, Rx, Rz) su qubit1, qubit2, qubit3
# --------------------------------------------------------------------------

@testset "3Q - SingleQubitGates" begin

    # Esempio: test su (X, Y, Z, H) e un range di θ per Rx, Rz
    # Per non esplodere, testiamoli su un solo bitstring di partenza (es. "000"),
    # o su 2-3 bitstring a scelta. Se vuoi testare TUTTI, puoi farlo.

    # 1) X, Y, Z, H su qubit1,2,3 partendo da "000"
    gates_1q = Dict(
        "X" => x_gate(),
        "Y" => y_gate(),
        "Z" => z_gate(),
        "H" => h_gate()
    )

    @testset "X, Y, Z, H on each qubit from |000>" begin
        for (gname, mat) in gates_1q
            for whichqubit in 1:3
                mps = create_product_mps("000")
                apply_single_qubit_gate!(mps, whichqubit, mat, 8)

                # costruiamo expected
                # apply (gname) to qubit 'whichqubit' in |000> => ...
                # Anziché derivare analiticamente, test "coerenza":
                #   per X: => se qubit=1 => |100>, qubit=2 => |010>, qubit=3 => |001>, ...
                #   ecc. Svolgiamo un if?

                expected = zeros(ComplexF64, 8)
                # index corrisponde a (a,b,c) => (qubit1, qubit2, qubit3), each in {0,1} => see for loop above.
                # con string "abc", dec index=4*a + 2*b + c + 1 se a,b,c in {0,1} interpretati in base2.
                # Piu semplice un if-elseif:

                if gname=="X"
                    # X|0> => |1>, X|1> => |0>
                    # => su "000", flip il qubit=whichqubit => bitstring ...
                    newbit = "100"  # se whichqubit=1
                    if whichqubit==2
                        newbit="010"
                    elseif whichqubit==3
                        newbit="001"
                    end
                    # scrivi la componente 1 in quell’indice
                    idx = bitstring_to_index3(newbit)
                    expected[idx] = 1
                elseif gname=="Y"
                    # Y|0> => i|1>, Y|1> => -i|0>
                    # qui su "000", => i*(|1> se whichqubit=1, etc.)
                    # => i su qubit=whichqubit => produce "100" con factor i?
                    newbit = "100"
                    if whichqubit==2
                        newbit="010"
                    elseif whichqubit==3
                        newbit="001"
                    end
                    idx = bitstring_to_index3(newbit)
                    expected[idx] = im
                elseif gname=="Z"
                    # Z|0> => +|0>, Z|1> => -|1>.
                    # su "000" => rimane "000" se the chosen qubit is in state 0 => factor +1
                    # => no change
                    # => expected = [1,0,0,0,0,0,0,0] => index=1
                    expected[1] = 1
                elseif gname=="H"
                    # H|0> => 1/sqrt(2)(|0> + |1>)
                    # => se qubit=1 => 1/sqrt(2)(|0,0,0> + |1,0,0>)
                    # => => index=1 + index=5 in base2 => (|000> + |100>)
                    val = 1/sqrt(2)
                    if whichqubit==1
                        idx0 = bitstring_to_index3("000")
                        idx1 = bitstring_to_index3("100")
                        expected[idx0] = val
                        expected[idx1] = val
                    elseif whichqubit==2
                        idx0 = bitstring_to_index3("000")
                        idx1 = bitstring_to_index3("010")
                        expected[idx0] = val
                        expected[idx1] = val
                    else
                        idx0 = bitstring_to_index3("000")
                        idx1 = bitstring_to_index3("001")
                        expected[idx0] = val
                        expected[idx1] = val
                    end
                end
                check_mps_state(mps, expected)
            end
        end
    end

    # 2) Rx, Rz su un set di thetas e su qubit1,2,3
    thetas = range(0, 2π, length=1000)
    for gatekind in ("Rx","Rz")
        @testset "$gatekind on qubit1,2,3" begin
            for whichqubit in 1:3
                for θ in thetas
                    mps = create_product_mps("000")
                    mat = (gatekind=="Rx") ? rx_gate(θ) : rz_gate(θ)
                    apply_single_qubit_gate!(mps, whichqubit, mat, 8)

                    # costruiamo expected
                    # come nel 2-qubit, ma i qubit rimanenti sono 0 => no change
                    # => se gatekind="Rx", su qubit=which => cos(θ/2)|0> - i sin(θ/2)|1>
                    # => qubit=1 => produce [cos(θ/2), 0,0,0, -i sin(θ/2), 0,0,0] ...
                    # e cosi via. Per brevità costruiamo vettore con una funzione di utility.

                    expected = build_expected_for_rxrz_3qubit(gatekind, whichqubit, θ, "000")
                    check_mps_state(mps, expected; atol=1e-7)
                end
            end
        end
    end

end # single-qubit gates

# Funzione di utilità per convertire bitstring 3 => index (1..8)
function bitstring_to_index3(bs::String)
    @assert length(bs)==3
    a = bs[1]=='1' ? 1 : 0
    b = bs[2]=='1' ? 1 : 0
    c = bs[3]=='1' ? 1 : 0
    # index in [1..8], 1 => a=b=c=0 => (0,0,0)
    return 4a + 2b + c + 1
end

# Esempio di costruttore per Rx/Rz "expected" su "000"
function build_expected_for_rxrz_3qubit(gkind, qubit, θ, initbits::String)
    # Per brevità, assumiamo initbits="000". Se volessi generico, dovresti estendere la logica
    vec = zeros(ComplexF64, 8)
    if gkind=="Rx"
        c = cos(θ/2)
        s = sin(θ/2)
        if qubit==1
            # => (c|0> - i s|1>) su qubit1, e qubit2=0, qubit3=0
            idx0 = bitstring_to_index3("000") # => c
            idx1 = bitstring_to_index3("100") # => -i s
            vec[idx0] = c
            vec[idx1] = -1im*s
        elseif qubit==2
            idx0 = bitstring_to_index3("000")
            idx1 = bitstring_to_index3("010")
            vec[idx0] = c
            vec[idx1] = -1im*s
        else
            idx0 = bitstring_to_index3("000")
            idx1 = bitstring_to_index3("001")
            vec[idx0] = c
            vec[idx1] = -1im*s
        end
    else
        # Rz => e^{-iθ/2}|0>, e^{+iθ/2}|1>
        # su "000" => qubit in state 0 => factor e^{-iθ/2} 
        # => global phase => [ e^{-iθ/2}, 0,0, ...] 
        # (stessa logica su qubit2=0 => no change, qubit3=0 => no change).
        # Manteniamo la phase => first element = e^{-iθ/2}
        # resto=0
        # NB: se qubit=1 => influisce su bit 1 => "000" rimane "000", con factor e^{-iθ/2}
        # se qubit=2 => => "000" => still "000", same factor e^{-iθ/2} if bit=0, etc.
        # In un test, e` sufficiente vedere se la contrazione produce la stessa wavefunction.
        # Mettiamo la phase in that position:
        vec[1] = exp(-1im*θ/2)
    end
    return vec
end

# --------------------------------------------------------------------------
# 4) Test Two-Qubit gates (SWAP, iSWAP, CNOT, CZ, ecc.)
#     su (1,2) e (2,3). (Se vuoi testare (1,3), potresti usare SWAP o aggiungere logica.)
# --------------------------------------------------------------------------

@testset "3Q - TwoQubitGates" begin

    # ----------------------------------------------------------------------
    # 1) Definiamo i gate a 2 qubit da testare e i pair (1,2) e (2,3)
    # ----------------------------------------------------------------------
    gates_2q = Dict(
        "SWAP" => swap_gate(),
        "iSWAP" => iswap_gate(),
        "CZ" => cz_gate(),
        "CNOT" => cnot_gate_clt()  # se definito in Gates.jl
    )
    pairs = [(1,2), (2,3)]

    # ----------------------------------------------------------------------
    # 3) Ciclo di test: per ogni gate2q, per ogni pair, per ogni bitstring
    # ----------------------------------------------------------------------

    for (gname, mat2q) in gates_2q
        @testset "Gate $gname su pairs e bitstring" begin
            for pair in pairs
                for bs in ALL_BITSTRINGS_3
                    # 1) Costruiamo lo stato MPS di partenza
                    mps = create_product_mps(bs)
                    # 2) Applichiamo gate2q su MPS
                    apply_two_qubit_gate!(mps, pair[1], mat2q, 8)

                    # 3) Costruiamo la vector 8D di partenza
                    initvec = build_3qubit_vector(bs)  # definisci sotto
                    # 4) Costruiamo l'operatore 8×8 e lo applichiamo
                    finalvec = apply_2qubit_op_global(mat2q, pair, initvec)

                    # 5) Confrontiamo con mps_to_vector_3qubit(mps)
                    check_mps_state(mps, finalvec)
                end
            end
        end
    end
end
    

# --------------------------------------------------------------------------
# 5) Test combinati (sequenze di gate singoli e doppi)
# --------------------------------------------------------------------------

@testset "3Q - Combined multi-step" begin

    # Funzioni per l'approccio globale (8D):
    # - build_3qubit_vector(bs) => Vector{ComplexF64} di dimensione 8 con 1 in pos corrispondente
    # - apply_single_qubit_op_global(mat2x2, pos, vec8) => newvec8
    # - apply_2qubit_op_global(mat4x4, (p1,p2), vec8) => newvec8
    #
    # E una utility per "aggiornare" la mappa (pos->qubit) quando facciamo SWAP(2,3), etc.

    """
        swap_positions!(mapping, p1, p2)

    Se 'mapping' è un array o dict che dice mapping[pos] = qubit,
    scambia i contenuti in p1 e p2.
    """
    function swap_positions!(mapping, p1, p2)
        temp = mapping[p1]
        mapping[p1] = mapping[p2]
        mapping[p2] = temp
    end

    # Caso A: spostiamo qubit3 vicino a qubit1 con SWAP(2,3), applichiamo iSWAP(1,2), Rz su pos2, X su pos3, e poi SWAP(2,3).
    @testset "Caso A bridging qubit3 -> pos2" begin
        for bs in ALL_BITSTRINGS_3
            # ----------------------------------------------------------------
            # 1) MPS approach
            # ----------------------------------------------------------------
            mpsA = create_product_mps(bs)

            # Step A1: SWAP(2,3) => scambia qubit2 e qubit3
            apply_two_qubit_gate!(mpsA, 2, swap_gate(), 8)

            # Step A2: iSWAP(1,2) => adesso in pos2 c'è qubit3
            apply_two_qubit_gate!(mpsA, 1, iswap_gate(), 8)

            # Step A3: Rz(pi/4) su pos2 (che adesso è qubit3)
            apply_single_qubit_gate!(mpsA, 2, rz_gate(pi/4), 8)

            # Step A4: X su pos3 (dove c'è qubit2)
            apply_single_qubit_gate!(mpsA, 3, x_gate(), 8)

            # Step A5: SWAP(2,3) => ripristina l'ordine (1,2,3)
            apply_two_qubit_gate!(mpsA, 2, swap_gate(), 8)

            # ----------------------------------------------------------------
            # 2) Global vector approach
            # ----------------------------------------------------------------
            # Prima, costruiamo la mappa: pos1-> qubit1, pos2-> qubit2, pos3-> qubit3
            # mapping[p] = qubit. Con 3 qubit, potremmo rappresentare con un array [1,2,3].
            mapping = [1,2,3]  # mapping[1]=1, mapping[2]=2, mapping[3]=3

            # init vec8
            vec = build_3qubit_vector(bs)

            # Step A1: SWAP(2,3) => scambiamo i qubit in pos2 e pos3
            #   global op => apply_2qubit_op_global(swap_gate(), (2,3), vec)
            #   e aggiorniamo mapping
            vec = apply_2qubit_op_global(swap_gate(), (2,3), vec)
            swap_positions!(mapping, 2, 3)

            # Step A2: iSWAP(1,2)
            #   phys. apply_2qubit_op_global(iSWAP, (1,2), vec)
            vec = apply_2qubit_op_global(iswap_gate(), (1,2), vec)

            # Step A3: Rz(pi/4) su pos2 (che e` mapping[2], v. qubit3?)
            #   devi fare apply_single_qubit_op_global(rz_gate(pi/4), pos=2, vec)
            vec = apply_single_qubit_op_global(rz_gate(pi/4), 2, vec)

            # Step A4: X su pos3 (che e` mapping[3], cioe` qubit2)
            vec = apply_single_qubit_op_global(x_gate(), 3, vec)

            # Step A5: di nuovo SWAP(2,3)
            vec = apply_2qubit_op_global(swap_gate(), (2,3), vec)
            swap_positions!(mapping, 2, 3)

            # Adesso 'mapping' dovrebbe tornare [1,2,3], e 'vec' e` lo stato finale globale.
            # Controlliamo l'MPS
            check_mps_state(mpsA, vec)
        end
    end

    # Caso B: spostiamo qubit1 a pos3 con SWAP(1,2), SWAP(2,3), poi applichiamo CZ(2,3), Y su pos2, X su pos3, e ritorniamo all'ordine
    @testset "Caso B bridging qubit1 -> pos3" begin
        for bs in ALL_BITSTRINGS_3
            # 1) MPS approach
            mpsB = create_product_mps(bs)
            apply_two_qubit_gate!(mpsB, 1, swap_gate(), 8)  # (1,2)
            apply_two_qubit_gate!(mpsB, 2, swap_gate(), 8)  # (2,3)
            # Ora l'ordine e` (2,3,1).
            apply_two_qubit_gate!(mpsB, 2, cz_gate(), 8)    # CZ(2,3)
            apply_single_qubit_gate!(mpsB, 2, y_gate(), 8)  # Y su pos2
            apply_single_qubit_gate!(mpsB, 3, x_gate(), 8)  # X su pos3
            # Riportiamo l'ordine: SWAP(2,3) -> (2,1,3), SWAP(1,2) -> (1,2,3)
            apply_two_qubit_gate!(mpsB, 2, swap_gate(), 8)
            apply_two_qubit_gate!(mpsB, 1, swap_gate(), 8)

            # 2) Global approach
            mapping = [1,2,3]
            vec = build_3qubit_vector(bs)

            # SWAP(1,2)
            vec = apply_2qubit_op_global(swap_gate(), (1,2), vec)
            swap_positions!(mapping, 1,2)  # adesso mapping=[2,1,3]
            # SWAP(2,3)
            vec = apply_2qubit_op_global(swap_gate(), (2,3), vec)
            swap_positions!(mapping, 2,3)  # mapping=[2,3,1]
            # CZ(2,3)
            vec = apply_2qubit_op_global(cz_gate(), (2,3), vec)
            # Y su pos2
            vec = apply_single_qubit_op_global(y_gate(), 2, vec)
            # X su pos3
            vec = apply_single_qubit_op_global(x_gate(), 3, vec)
            # poi SWAP(2,3), SWAP(1,2) per tornare [1,2,3]
            vec = apply_2qubit_op_global(swap_gate(), (2,3), vec)
            swap_positions!(mapping, 2,3)  # => [2,1,3]
            vec = apply_2qubit_op_global(swap_gate(), (1,2), vec)
            swap_positions!(mapping, 1,2)  # => [1,2,3]

            # Controllo finale
            check_mps_state(mpsB, vec)
        end
    end

    # Caso C: sequenza più lunga (esempio) partendo da "111", spostiamo qubit3 in pos2, CNOT(1->2), Rz e Rx su pos2, X su pos3, e riportiamo l'ordine
    @testset "Caso C bridging e sequenza lunga" begin
        for bs in ALL_BITSTRINGS_3
            # 1) MPS
            mpsC = create_product_mps(bs)
            # SWAP(2,3)
            apply_two_qubit_gate!(mpsC, 2, swap_gate(), 8)
            # CNOT(1->2)
            apply_two_qubit_gate!(mpsC, 1, cnot_gate_clt(), 8)
            # Rz(pi/6) su pos2
            apply_single_qubit_gate!(mpsC, 2, rz_gate(pi/6), 8)
            # Rx(pi/2) su pos2
            apply_single_qubit_gate!(mpsC, 2, rx_gate(pi/2), 8)
            # X su pos3
            apply_single_qubit_gate!(mpsC, 3, x_gate(), 8)
            # Ri-swap(2,3)
            apply_two_qubit_gate!(mpsC, 2, swap_gate(), 8)

            # 2) Globale
            mapping = [1,2,3]
            vec = build_3qubit_vector(bs)
            # swap(2,3)
            vec = apply_2qubit_op_global(swap_gate(), (2,3), vec)
            swap_positions!(mapping,2,3)
            # cnot(1->2)
            vec = apply_2qubit_op_global(cnot_gate_clt(), (1,2), vec)
            # Rz(pi/6) su pos2
            vec = apply_single_qubit_op_global(rz_gate(pi/6), 2, vec)
            # Rx(pi/2) su pos2
            vec = apply_single_qubit_op_global(rx_gate(pi/2), 2, vec)
            # X su pos3
            vec = apply_single_qubit_op_global(x_gate(), 3, vec)
            # swap(2,3)
            vec = apply_2qubit_op_global(swap_gate(), (2,3), vec)
            swap_positions!(mapping,2,3)

            # check
            check_mps_state(mpsC, vec)
        end
    end

end
    
