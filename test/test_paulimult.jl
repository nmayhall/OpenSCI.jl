using PauliOperators
using OpenSCI
using Test
using LinearAlgebra
using Printf
using Random
using Plots

@testset "Paulimult" begin
# function test1()
    Random.seed!(1234)
    N = 3

    L = Lindbladian(N)
    add_hamiltonian!(L, OpenSCI.heisenberg_1D(N, 1.1, 1.2, 1.3))
    add_channel_dephasing!(L, .1)
    add_channel_depolarizing!(L, .1)
    add_channel_amplitude_damping!(L, .1)
    display(L)

    o = rand(Pauli{N})
    display(o)

    println(" Test Lindbladian times PauliSum")
    err = Matrix(L)*vec(Matrix(o)) - vec(Matrix(L*o))
    @test norm(err) < 1e-14
    

    println(" Test Lindbladian adjoint times PauliSum")
    err = Matrix(L)'*vec(Matrix(o)) - vec(Matrix(L'*o))
    @test norm(err) < 1e-14
    return 
end

# test1()