module States

using LinearAlgebra

export ket0, ket1, ketp, ketm,
       ket2, ket3,  # stati di 2 qubit base
       kron2,       # per comporre stati a 2 qubit generici
       normdiff     # funzione di confronto

"""
    ket0() -> Vector{ComplexF64}

Ritorna il vettore |0> in dimensione 2.
"""
function ket0()::Vector{ComplexF64}
    return [1+0im, 0+0im]
end

"""
    ket1() -> Vector{ComplexF64}

Ritorna il vettore |1> in dimensione 2.
"""
function ket1()::Vector{ComplexF64}
    return [0+0im, 1+0im]
end

"""
    ketp() -> Vector{ComplexF64}

Ritorna |+> = 1/sqrt(2) [|0> + |1>].
"""
function ketp()::Vector{ComplexF64}
    return (1/sqrt(2)) * [1+0im, 1+0im]
end

"""
    ketm() -> Vector{ComplexF64}

Ritorna |-> = 1/sqrt(2) [|0> - |1>].
"""
function ketm()::Vector{ComplexF64}
    return (1/sqrt(2)) * [1+0im, -1+0im]
end

"""
    kron2(u, v) -> Vector{ComplexF64}

Restituisce il prodotto di Kronecker tra i vettori u e v (entrambi 2D).
Rappresenta lo stato a 2 qubit |uv>.
"""
function kron2(u::Vector{<:ComplexF64}, v::Vector{<:ComplexF64})
    return kron(u, v)  # builtin in LinearAlgebra
end

"""
    normdiff(x, y) -> Real

Restituisce la norma Euclidea di (x - y), utile per testare la vicinanza di stati.
"""
function normdiff(x, y)
    return norm(x .- y)
end

end # module
