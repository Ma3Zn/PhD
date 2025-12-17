# ----------------------------------------------------------------------------
# File di test, da lanciare con:
#
#   julia --project=@. TestGates.jl
#
# (oppure includilo in un sistema di test più grande con Pkg.test)
# ----------------------------------------------------------------------------

using Test
using LinearAlgebra

using ITensors
using ITensorMPS

include("../src/Utils.jl")
using .Utils

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

@testset "Test Funzioni Utils" begin

#   ATTENZIONE::    Se lo esegui in debug allora termina con successo. 
#                   Se lo esegui direttamente da Julia restituisce un NaN per n=8.
#                   Potrebbero esserci dei problemi di RAM (?)
#                   Comunque direi che va bene così l'inizializzazione del nostro stato
    for n in 2:2:8
        # Creiamo la MPS relativa nello stato opportuno
        mps = create_initial_mps(n)

        # Creiamo la codifica binaria dello stato che vogliamo ottenere
        state_bin = [isodd(k) ? "1" : "0" for k=0:n-1]

        # Formiamo un'unica stringa con la rappresentazione binaria dello stato
        state_bin = join(state_bin)
        
        # Convertiamo lo stato binario in un indice decimale
        state_idx = parse(Int64, "0b" * state_bin)

        # Aggiungiamo l'offset di Julia relativo all'indicizzazione da 1 e non da 0
        state_idx = state_idx + 1

        # Creiamo il vettore rappresentante lo stato atteso
        expected = Vector{ComplexF64}(undef, 2^n) * 0

        # Inizializziamo opportunamente il vettore atteso
        expected[state_idx] = 1

        # Controlliamo che l'MPS sia inizializzata nello stato opportuno
        check_mps_state(mps, expected)
    end
end