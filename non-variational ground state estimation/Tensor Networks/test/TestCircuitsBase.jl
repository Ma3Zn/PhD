using Test
using ITensors
using ITensorMPS

# Includi i file necessari:
include("../src/Circuit.jl")  # definisce apply_single_qubit_gate!, apply_two_qubit_gate!
using .Circuit

include("../src/Gates.jl")    # definisce x_gate, h_gate, cnot_gate, swap_gate, ecc.
using .Gates

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


###############################################################################
# Test Suite
###############################################################################

@testset "Test apply_single_qubit_gate! e apply_two_qubit_gate!" begin

#   ATTENZIONE::    l'ordine dei qubit è da sinistra a destra 
#                             --> |Q1,Q2,...,Qn>

    # Impostiamo un chi ragionevole
    chi = 4

    # -------------------------------------------------------------------------
    # CASO 2 SITI
    # -------------------------------------------------------------------------
    @testset "2 SITI" begin
        # Creo MPS = |00>
        mps_2 = create_product_mps("00")

        # 1) Applico X al sito 1 => stato atteso = |10>
        apply_single_qubit_gate!(mps_2, 1, x_gate(), chi)
        expected = [0, 0, 1, 0]  # in ordine |00>, |01>, |10>, |11> --- |Q1,Q2>
        check_mps_state(mps_2, expected)

        # 2) Applico X al sito 2 => |11>
        apply_single_qubit_gate!(mps_2, 2, x_gate(), chi)
        expected = [0, 0, 0, 1]
        check_mps_state(mps_2, expected)

        # 3) Riporto lo stato a |00> per test successivi
        #    (Applico X su sito 1 e 2)
        apply_single_qubit_gate!(mps_2, 1, x_gate(), chi)
        apply_single_qubit_gate!(mps_2, 2, x_gate(), chi)
        # Stato atteso di nuovo [1,0,0,0]
        check_mps_state(mps_2, [1, 0, 0, 0])

        # 4) Applico H sul sito 1 => (|0> + |1>)/√2 tensore |0>
        apply_single_qubit_gate!(mps_2, 1, h_gate(), chi)
        # Ora lo stato = 1/sqrt(2)[ |00> + |10> ]
        # => in vettore = [1/sqrt(2), 0, 1/sqrt(2), 0]
        val = 1/sqrt(2)
        expected = [val, 0, val, 0]
        check_mps_state(mps_2, expected)

        # 5) Applico SWAP(1,2)
        apply_two_qubit_gate!(mps_2, 1, swap_gate(), chi)
        # scambia i 2 qubit => (1/sqrt(2)) [|00> + |01>] 
        expected = [val, val, 0, 0]
        check_mps_state(mps_2, expected)

        # 6) Applico SWAP(1,2)
        apply_two_qubit_gate!(mps_2, 1, swap_gate(), chi)
        # scambia i 2 qubit => (1/sqrt(2)) [|00> + |10>] 
        expected = [val, 0, val, 0]
        check_mps_state(mps_2, expected)

        # 7) Applico SWAP(1,2)
        apply_two_qubit_gate!(mps_2, 1, swap_gate(), chi)
        # scambia i 2 qubit => (1/sqrt(2)) [|00> + |01>] 
        expected = [val, val, 0, 0]
        check_mps_state(mps_2, expected)

        # 8) Eseguo una decomposizione dello SWAP in 3 CNOT
        apply_two_qubit_gate!(mps_2, 1, cnot_gate_clt(), chi)
        apply_two_qubit_gate!(mps_2, 1, cnot_gate_cgt(), chi)
        apply_two_qubit_gate!(mps_2, 1, cnot_gate_clt(), chi)
        # scambia i 2 qubit => (1/sqrt(2)) [|00> + |10>] 
        expected = [val, 0, val, 0]
        check_mps_state(mps_2, expected)

        # 9) Applico CNOT(1->2)
        apply_two_qubit_gate!(mps_2, 1, cnot_gate_clt(), chi)
        # (1/sqrt(2)) [|00> + |11>]
        expected = [val, 0, 0, val]
        check_mps_state(mps_2, expected)
    end

    # -------------------------------------------------------------------------
    # CASO 3 SITI
    # -------------------------------------------------------------------------
    @testset "3 SITI" begin
        # Creo MPS = |000>
        mps_3 = create_product_mps("000")

        # 1) Applico X su sito 3 => |001>
        apply_single_qubit_gate!(mps_3, 3, x_gate(), chi)
        # Vettore atteso: dimensione 2^3=8
        # Ordinamento: |000>,|001>,|010>,|011>,|100>,|101>,|110>,|111>
        expected = [0,1,0,0,0,0,0,0]
        check_mps_state(mps_3, expected)

        # 2) Applico X su sito 1 => => |100>
        apply_single_qubit_gate!(mps_3, 1, x_gate(), chi)
        # adesso stiamo modificando dal vettore |001>
        # => X su qubit 1 => sposta |001> in |101>, giusto?
        # Facciamo attenzione all’ordinamento qubit(1), qubit(2), qubit(3).
        # Se qubit(1) è la posizione "più a sinistra", allora |001> => (qubit1=0, qubit2=0, qubit3=1).
        # Applicare X su qubit1 => (qubit1=1, qubit2=0, qubit3=1) => |101>.
        # In vettore (|000>, |001>, |010>, |011>, |100>, |101>, |110>, |111>),
        # |101> corrisponde a 6° entry => [0,0,0,0,0,1,0,0]
        expected = [0,0,0,0,0,1,0,0]
        check_mps_state(mps_3, expected)

        # 3) Applico CNOT(2->3) (sites 2,3)
        #    Siti adiacenti i=2 => qubit2, i+1=3 => qubit3
        apply_two_qubit_gate!(mps_3, 2, cnot_gate_clt(), chi)
        # Lo stato era |101>, qubit2=0 => target non flippa => resta |101>
        expected = [0,0,0,0,0,1,0,0]
        check_mps_state(mps_3, expected)

        # 4) Applico X su qubit2 => lo stato diventa |111>
        apply_single_qubit_gate!(mps_3, 2, x_gate(), chi)
        # |101> => qubit2=0 -> 1 => |111>
        # => [0,0,0,0,0,0,0,1]
        expected = [0,0,0,0,0,0,0,1]
        check_mps_state(mps_3, expected)

        # 5) Applico SWAP(1,2)
        #    Scambia i qubit 1 e 2 => |111> => |111>, in base standard non cambia nulla
        apply_two_qubit_gate!(mps_3, 1, swap_gate(), chi)
        check_mps_state(mps_3, expected)

        # TODO:: Estendere abbondantemente
    end

    # -------------------------------------------------------------------------
    # CASO 4 SITI
    # -------------------------------------------------------------------------
    @testset "4 SITI" begin
        # MPS = |0000>
        mps_4 = create_product_mps("0000")

        # 1) Applico H su qubit 1 => (|0> + |1>)/√2 su qubit1, e |000> sugli altri
        apply_single_qubit_gate!(mps_4, 1, h_gate(), chi)
        # lo stato = 1/sqrt(2)[|0000> + |1000>]
        # => in vettore dimensione 16, index 1 => 0..15
        #  |0000> => indice 0 => coeff = 1/sqrt(2)
        #  |1000> => indice 8 => coeff = 1/sqrt(2)
        #  => [1/sqrt(2),0,0,0,0,0,0,0, 1/sqrt(2),0,0,0,0,0,0,0]
        val = 1/sqrt(2)
        expected = zeros(ComplexF64, 16)
        expected[1]   = val  # 1-based indexing => expected[1] = |0000>
        expected[9]   = val  # => |1000> corrisponde a index=9
        check_mps_state(mps_4, expected)

        # 2) Applico CNOT(1->2)
        #    => se qubit1=1 => flip qubit2
        apply_two_qubit_gate!(mps_4, 1, cnot_gate_clt(), chi)
        # Quella parte era |1000> => qubit1=1, qubit2=0 => qubit2=1 => |1100>
        # Quindi lo stato diventa 1/sqrt(2)[|0000> + |1100>]
        # => index(12) = |1100> => 1-based: offset=12+1=13
        fill!(expected, 0)
        expected[1]  = val   # |0000>
        expected[13] = val   # |1100> => binario 1100 => dec 12 => pos 13
        check_mps_state(mps_4, expected)

        # TODO:: Estendere ancora piu' abbondantemente
    end

end
