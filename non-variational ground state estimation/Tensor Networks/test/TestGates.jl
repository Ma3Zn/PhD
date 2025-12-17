# ----------------------------------------------------------------------------
# File di test, da lanciare con:
#
#   julia --project=@. TestGates.jl
#
# (oppure includilo in un sistema di test più grande con Pkg.test)
# ----------------------------------------------------------------------------

using Test
using LinearAlgebra

include("../src/Gates.jl")
using .Gates

include("../src/States.jl")
using .States

@testset "Test Gates Matriciali" begin

    # Identità 2×2 e 4×4 per i test
    I2 = Matrix{Float64}(I, 2, 2)
    I4 = Matrix{Float64}(I, 4, 4)

    # 1) X, Y, Z
    @test x_gate()*x_gate() ≈ I2
    @test y_gate()*y_gate() ≈ I2
    @test z_gate()*z_gate() ≈ I2

    # 2) Hadamard
    # Verifichiamo che H^2 = I (fino a possibili fattori numerici)
    H = h_gate()
    @test H*H ≈ I2

    # 3) CNOT
    # cnot^2 = I4
    cnot = cnot_gate()
    @test cnot*cnot ≈ I4

    # 4) SWAP
    # swap^2 = I4
    swap = swap_gate()
    @test swap*swap ≈ I4

    # 5) Rotazioni X e Z
    #   Verifichiamo che Rx(theta)*Rx(-theta) = I
    #   e Rz(theta)*Rz(-theta) = I, con una certa tolleranza
    θ = 0.3
    Rx = rx_gate(θ)
    Rx_inv = rx_gate(-θ)
    @test Rx*Rx_inv ≈ Matrix{ComplexF64}(I, 2, 2)

    Rz = rz_gate(θ)
    Rz_inv = rz_gate(-θ)
    @test Rz*Rz_inv ≈ Matrix{ComplexF64}(I, 2, 2)

#   ---------------------------
#   Sezione 2: test azione su stati a singolo qubit
#   ---------------------------
#   X|0> => |1>
    v0 = ket0()
    v1 = ket1()
    @test normdiff(x_gate()*v0, v1) < 1e-10
    @test normdiff(x_gate()*v1, v0) < 1e-10

#   Y|0> => i|1>, Y|1> => -i|0>
    @test normdiff(y_gate()*v0, im*v1) < 1e-10
    @test normdiff(y_gate()*v1, -im*v0) < 1e-10

#   Z|0> => |0>, Z|1> => -|1>
    @test normdiff(z_gate()*v0, v0) < 1e-10
    @test normdiff(z_gate()*v1, -v1) < 1e-10

#   H|0> => |+>, H|+> => |0>, e simili
    vp = ketp()
    @test normdiff(h_gate()*v0, vp) < 1e-10
    @test normdiff(h_gate()*vp, v0) < 1e-10

#   ---------------------------
#   Sezione 3: test azione su stati a 2 qubit
#   ---------------------------
#   3.1) CNOT
#    CNOT|00> => |00>, CNOT|01> => |01>, CNOT|10> => |11>, CNOT|11> => |10>
    ket00 = kron2(v0, v0)
    ket01 = kron2(v0, v1)
    ket10 = kron2(v1, v0)
    ket11 = kron2(v1, v1)

    @test normdiff(cnot*ket00, ket00) < 1e-10
    @test normdiff(cnot*ket01, ket01) < 1e-10
    @test normdiff(cnot*ket10, ket11) < 1e-10
    @test normdiff(cnot*ket11, ket10) < 1e-10

#   3.2) SWAP
#    SWAP|01> => |10>, SWAP|10> => |01>, etc.
    @test normdiff(swap*ket01, ket10) < 1e-10
    @test normdiff(swap*ket10, ket01) < 1e-10
    @test normdiff(swap*ket00, ket00) < 1e-10
    @test normdiff(swap*ket11, ket11) < 1e-10

#   3.3) CNOT su |+0> => (1/sqrt(2)) [|00> + |10>] => (1/sqrt(2)) [|00> + |11>]
    ketp0 = kron2(vp, v0)
    expected = (1/sqrt(2)) * (ket00 .+ ket11)
    @test normdiff(cnot*ketp0, expected) < 1e-10

#   3.4) SWAP su |+0> => |0+> (cioè scambia i due qubit)
    ket0p = kron2(v0, vp)
    @test normdiff(swap*ketp0, ket0p) < 1e-10

end
