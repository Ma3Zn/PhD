module Utils

# ======================================================
#  Utils.jl
#
#  Contiene parametri e funzioni di supporto per
#  la simulazione (periodicità, costanti d'interazione,
#  funzioni di drive lineare e single-particle).
# ======================================================

export  closed_ring, J, lam, lam_sp, drive, drive_sp, 
        create_initial_mps, diag_sqrt!, qubit_siteind

using Base: pi
using LinearAlgebra
using ITensors
using ITensorMPS

"""
    closed_ring :: Bool

Se `true`, indica connettività periodica (anello).
"""
const closed_ring = true

"""
    J :: NTuple{3,Float64}

Costanti d'interazione [Jx, Jy, Jz] in unità di J.
"""
const J = [1.0, 1.0, 1.0]

"""
    lam :: Float64

Parametro che regola la velocità del processo di annealing.
"""
const lam = 0.1

"""
    lam :: Float64

Parametro che regola la velocità di spegnimento del termine a singolo corpo.
"""
const lam_sp = lam

"""
    drive(param::Real, t::Real) -> Real

Drive lineare: ris = param * t

Drive sinusoidale: ris = sin((param*t)*pi/2)^0.9
"""
function drive(param::Real, t::Real)::Float64
    # Linear
    # ris = param * t

    # Esempio alternativo (sinusoidale):
    # ris = (sin((param*t)*pi/2))^0.9

    return sin((param*t)*pi/2)^0.9
end

"""
    drive_sp(param::Real, n::Real, t::Real) -> Float64

Drive a singolo corpo, formula "potenza":
  se (param * t)^n < 1 allora ris = 1 - (param * t)^n
  altrimenti ris = 0

formula "tangente":
    tan(1) - tan((param * t)^n)
"""
function drive_sp(param::Real, n::Real, t::Real)::Float64
    # val = (param * t)^n
    # if val < 1
    #     return 1 - val
    # else
    #     return 0.0
    # end

    return tan(1) - tan((param * t)^n)
end

function create_initial_mps(N::Int)
"""
        create_initial_mps(N) -> MPS
    
    Crea un MPS di N qubit nello stato |000...0>.
"""

#   Crea un "site set" di tipo S=1/2, dimensione = 2 per ciascun sito
    sites = siteinds("S=1/2", N)

#   Creiamo il ground state dell'Hamiltoniana iniziale
#   Up --> |0>, Dn --> |1>
    state = [isodd(n) ? "Dn" : "Up" for n=0:N-1]

#   Costruiamo un MPS nello stato desiderato
    return MPS(sites, state)

end

"""
        diag_sqrt -> ITensor

    Dato un MPS diagonale[!] ne calcola la radice quadrata
"""
function diag_sqrt!(S::ITensor)
    # Supponiamo che S abbia esattamente 2 indici link, 
    # ad es. (l, l') se è il tensore diagonale dall'SVD.
    indsS = inds(S)
    @assert length(indsS) == 2 "Si assume che S abbia 2 indici (left, right)."
    iL, iR = indsS[1], indsS[2]
    
    # Iteriamo su tutti i possibili valori dell'indice iL (dimension = dim(iL))
    # e prendiamo la componente diagonale (iL=>n, iR=>n).
    for n in 1:dim(iL)
        val = S[iL => n, iR => n]
        S[iL => n, iR => n] = sqrt(val)
    end
    return S
end

"""
    siteind(mps, i) -> Index

Ritorna l'Index fisico del sito i, cercando fra gli indici di mps[i]
quello che possiede i tag "Qubit", "Site", e "n=i".

ATTENZIONE:: inutile se usi il tag S=1/2 per i tuoi siti, l'interfaccia
             siteind fornita da ITensors funziona già perfettamente
"""
function qubit_siteind(mps::MPS, i::Int)
    Ai = mps[i]
    allI = inds(Ai)
    # cerca l'indice che ha tag "Qubit,Site"
    found = filter(idx -> hastags(idx, "S=1/2"), allI)

    if length(found) == 1
        return found[1]
    else
        error("siteind(mps, $i): nessun indice con tag Qubit in mps[$i]!")
    end
end

end # module
