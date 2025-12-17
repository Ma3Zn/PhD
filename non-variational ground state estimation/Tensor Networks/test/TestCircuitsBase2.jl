using Test
using ITensors
using ITensorMPS

# Se hai i moduli Circuit, Gates, States, etc., includi:
include("../src/Circuit.jl")   # se definisce apply_single_qubit_gate!, apply_two_qubit_gate!
using .Circuit

include("../src/Gates.jl")     # definizioni X, Y, Z, H, rx_gate, rz_gate, swap_gate, iswap_gate, etc.
using .Gates

###############################################################################
# Sezione 1: utility specifiche per test su 2 qubit
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
# Sezione 1: test su singolo qubit (H, Y, Z, Rx, Rz)
###############################################################################

@testset "2 qubit - SingleQubitGates" begin

    # |00>, |01>, |1,0>, |1,1>

    # 1) H test
    @testset "Apply H on qubit1, qubit2" begin
        # MPS = |00>
        for whichqubit in (1,2)
            mps = create_product_mps("00")
            # apply H on qubit=whichqubit
            apply_single_qubit_gate!(mps, whichqubit, h_gate(), 4)  # chi=4
            # expected = H|0> tensore |0> (se whichqubit=1) oppure |0> tensore H|0> (if 2)
            # nel primo caso => 1/sqrt(2)[|0,0> + |1,0>]
            # => vector = [1/sqrt(2), 0, 1/sqrt(2), 0]
            # se qubit2 => 1/sqrt(2)[|0,0> + |0,1>]
            # => vector = [1/sqrt(2), 1/sqrt(2), 0, 0]
            vec = zeros(ComplexF64, 4)
            val = 1/sqrt(2)
            if whichqubit == 1
                vec[1] = val
                vec[3] = val
            else
                vec[1] = val
                vec[2] = val
            end
            check_mps_state(mps, vec)
        end
    end

    # 2) Y, Z
    @testset "Apply Y, Z on qubit1" begin
        # MPS = |00>
        mps = create_product_mps("00")
        apply_single_qubit_gate!(mps, 1, y_gate(), 4)
        # Y|0> => i|1>, qubit1 => i|1>, qubit2 => |0>
        # => final = i|10> => vector = [0,0,i,0]
        # check
        vec = [0+0im, 0+0im, im, 0+0im]
        check_mps_state(mps, vec)

        # Applico Z sul qubit2 => da i|10> => i|10> (since qubit2=0 => Z|0> = +1|0>)
        apply_single_qubit_gate!(mps, 2, z_gate(), 4)
        check_mps_state(mps, vec)
    end

    # 3) Rx, Rz su una gamma di theta
    @testset "Rx,Rz param test" begin
        using LinearAlgebra: norm
        thetas = range(0, 2π, length=1000)  # 7 valori
        for whichqubit in (1,2)
            for gatekind in ("Rx", "Rz")
                for θ in thetas
                    # MPS = |00>
                    mps = create_product_mps("00")
                    # costruiamo la matrice
                    mat = gatekind=="Rx" ? rx_gate(θ) : rz_gate(θ)

                    # Apply
                    apply_single_qubit_gate!(mps, whichqubit, mat, 4)

                    # Confronto con formula analitica: 
                    # if gatekind=="Rx" => Rx(θ)|0> = cos(θ/2)|0> - i sin(θ/2)|1>
                    # if gatekind=="Rz" => Rz(θ)|0> = e^{-iθ/2} (cos(θ/2)|0> - sin(θ/2)|0> ???)
                    # Actually Rz(θ)|0> => exp(-iθ/2)|0>, Rz(θ)|1> => exp(+iθ/2)|1>.
                    # Su |0> => e^{-iθ/2}|0>. 
                    # Poi c'è qubit2 => |0> invariato.

                    # costruiamo expected vector
                    vec = zeros(ComplexF64, 4)
                    if whichqubit == 1
                        # qubit1 subisce gate, qubit2 rimane |0>
                        # => dimension: 
                        #   if Rx(θ):
                        #     cos(θ/2)|0> - i sin(θ/2)|1>  tensore  |0>
                        # => [ cos(θ/2), 0, -i sin(θ/2), 0 ]
                        #   if Rz(θ):
                        #     e^{-iθ/2}|0>  => => [ e^{-iθ/2}, 0, 0, 0 ]
                        #  (il gate su |1> => e^{+iθ/2}|1>, ma qubit1=|0> di partenza.)
                        if gatekind=="Rx"
                            c = cos(θ/2)
                            s = sin(θ/2)
                            vec[1] = c + 0im
                            vec[3] = -1im * s
                        else
                            # Rz(θ) on |0> => e^{-iθ/2}|0>
                            # vector => [ e^{-iθ/2}, 0, 0, 0 ]
                            # = e^{-iθ/2}|0> => global phase might matter. 
                            # spostiamo la global phase se vogliamo "real" (non definita)
                            # Teniamola come e^{-iθ/2}.
                            vec[1] = exp(-1im*θ/2)
                        end
                    else
                        # qubit2 subisce gate, qubit1=|0>.
                        # =>  qubit1 => |0>, qubit2 => gate|0>.
                        # => dimension:
                        #   if Rx: cos(θ/2)|00> - i sin(θ/2)|01>
                        # => [ cos(θ/2), -i sin(θ/2), 0, 0 ]
                        if gatekind=="Rx"
                            c = cos(θ/2)
                            s = sin(θ/2)
                            vec[1] = c
                            vec[2] = -1im * s
                        else
                            # Rz(θ) on qubit2=|0> => e^{-iθ/2}|0>, 
                            # => [ e^{-iθ/2}, 0, 0, 0 ]
                            vec[1] = exp(-1im*θ/2)
                        end
                    end

                    check_mps_state(mps, vec; atol=1e-7)
                end
            end
        end
    end

end # testset single-qubit gates

###############################################################################
# Sezione 2: test gate a due qubit (nuovi)
###############################################################################

@testset "2 qubit - TwoQubitGates" begin

    # 1) SWAP su |01> => |10>, etc.
    @testset "SWAP gate" begin
        for bitstring in ("00","01","10","11")
            mps = create_product_mps(bitstring)
            apply_two_qubit_gate!(mps, 1, swap_gate(), 4)
            # expected => scambia qubit1 e qubit2
            # "01" => "10", "10" => "01", "00"=> "00", "11"=>"11"
            swapped = if bitstring == "01"
                "10"
            elseif bitstring == "10"
                "01"
            else
                bitstring
            end
            # costruiamo expected vector
            # mps corrispondente a swapped
            mps_expected = create_product_mps(swapped)
            vec_expected = mps_to_vector(mps_expected)
            check_mps_state(mps, vec_expected)
        end
    end

    # 2) iSWAP
    @testset "iSWAP gate" begin
        for bitstring in ("00","01","10","11")
            mps = create_product_mps(bitstring)
            apply_two_qubit_gate!(mps, 1, iswap_gate(), 4)
            # iSWAP => 
            #   |00> -> |00>,
            #   |11> -> |11>,
            #   |01> -> i|10>,
            #   |10> -> i|01>.
            # costruiamo expected
            vec = zeros(ComplexF64, 4)
            if bitstring == "00"
                vec[1] = 1
            elseif bitstring == "01"
                # i|10> => i*(0,0,1,0)
                vec[3] = im
            elseif bitstring == "10"
                vec[2] = im
            else # "11"
                vec[4] = 1
            end
            check_mps_state(mps, vec)
        end
    end

    # 3) CZ
    @testset "CZ gate" begin
        for bitstring in ("00","01","10","11")
            mps = create_product_mps(bitstring)
            apply_two_qubit_gate!(mps, 1, cz_gate(), 4)
            # CZ => |00>->|00>, |01>->|01>, |10>->|10>, |11>-> -|11>
            vec = zeros(ComplexF64,4)
            if bitstring == "00"
                vec[1] = 1
            elseif bitstring == "01"
                vec[2] = 1
            elseif bitstring == "10"
                vec[3] = 1
            else
                vec[4] = -1
            end
            check_mps_state(mps, vec)
        end
    end

end # testset two-qubit gates
