module Gates

# Qui definiamo i gate in forma di Array Julia.
# Eventualmente puoi trasformarli in ITensor quando li applichi all'MPS.

export x_gate, y_gate, z_gate, h_gate, rx_gate, rz_gate, cnot_gate_clt, cnot_gate_cgt, swap_gate, iswap_gate, cz_gate

using LinearAlgebra

"""
    x_gate() :: Matrix{Float64}

Matrice 2×2 per il gate di Pauli-X.
"""
function x_gate()::Matrix{Float64}
    return [0 1;
            1 0]
end

"""
    y_gate() :: Matrix{ComplexF64}

Matrice 2×2 per il gate di Pauli-Y.
"""
function y_gate()::Matrix{ComplexF64}
    return [0 -im;
            im  0]
end

"""
    z_gate() :: Matrix{Float64}

Matrice 2×2 per il gate di Pauli-Z.
"""
function z_gate()::Matrix{Float64}
    return [1  0;
            0 -1]
end

"""
    h_gate() :: Matrix{Float64}

Matrice 2×2 per la Hadamard, con fattore 1/√2.
"""
function h_gate()::Matrix{Float64}
    return (1/sqrt(2)) * [1  1;
                          1 -1]
end

"""
    rx_gate(theta::Real) :: Matrix{ComplexF64}

Matrice 2×2 per la rotazione attorno all'asse X di un angolo theta.
Rx(θ) = cos(θ/2)*I - i*sin(θ/2)*X
"""
function rx_gate(theta::Real)::Matrix{ComplexF64}
    c = cos(theta/2)
    s = sin(theta/2)
    #  [ c    -i s ]
    #  [-i s   c   ]
    return [ c   -1im*s;
            -1im*s   c  ]
end

"""
    rz_gate(theta::Real) :: Matrix{ComplexF64}

Matrice 2×2 per la rotazione attorno all'asse Z di un angolo theta.
Rz(θ) = cos(θ/2)*I - i*sin(θ/2)*Z
"""
function rz_gate(theta::Real)::Matrix{ComplexF64}
    c = cos(theta/2)
    s = sin(theta/2)
    #  [ c - i s      0       ]
    #  [     0    c + i s     ]
    return [ c - 1im*s     0;
             0        c + 1im*s ]
end

"""
    cnot_gate() :: Matrix{Float64}

Matrice 4×4 per la CNOT (ctrl=qubit 1, targ=qubit 2).
Con convenzione: |00>, |01>, |10>, |11>.

ATTENZIONE::    L'ordinamento dei qubit è da sinistra a destra
                        ---> |Q1,Q2,...,Qn]
"""
function cnot_gate_clt()::Matrix{Float64}
#   Qc < Qt => |Qc,Qt>
    return [ 1  0  0  0;
             0  1  0  0;
             0  0  0  1;
             0  0  1  0 ]
end

"""
    cnot_gate() :: Matrix{Float64}

Matrice 4×4 per la CNOT (ctrl=qubit 2, targ=qubit 1).
Con convenzione: |00>, |01>, |10>, |11>.

ATTENZIONE::    L'ordinamento dei qubit è da sinistra a destra
                        ---> |Q1,Q2,...,Qn]
"""
function cnot_gate_cgt()::Matrix{Float64}
#   Qc > Qt => |Qt,Qc>
    return [ 1  0  0  0;
             0  0  0  1;
             0  0  1  0;
             0  1  0  0 ]
end

"""
    swap_gate() :: Matrix{Float64}

Matrice 4×4 per lo swap (ctrl=qubit 1, targ=qubit 2).
Con convenzione: |00>, |01>, |10>, |11>.
"""
function swap_gate()::Matrix{Float64}
    return [ 1  0  0  0;
             0  0  1  0;
             0  1  0  0;
             0  0  0  1 ]
end

"""
    iswap_gate() :: Matrix{Float64}

Matrice 4×4 per lo iswap (ctrl=qubit 1, targ=qubit 2).
Con convenzione: |00>, |01>, |10>, |11>.
"""
function iswap_gate()::Matrix{ComplexF64}
    # iSWAP scambia |01> con |10> e aggiunge fattore i
    # la matrice standard:
    #   [1  0  0  0
    #    0  0  i  0
    #    0  i  0  0
    #    0  0  0  1]
    M = zeros(ComplexF64, 4,4)
    M[1,1] = 1
    M[4,4] = 1
    M[2,3] = im
    M[3,2] = im
    return M
end

"""
    CZ_gate() :: Matrix{Float64}

Matrice 4×4 per il CZ (ctrl=qubit 1, targ=qubit 2).
Con convenzione: |00>, |01>, |10>, |11>.
"""
function cz_gate()::Matrix{Float64}
    # |00> => 1, |01> => 1, |10> => 1, |11> => -1
    # diagonale(1,1,1,-1)
    return [1 0 0 0;
            0 1 0 0;
            0 0 1 0;
            0 0 0 -1]
end

end # module
