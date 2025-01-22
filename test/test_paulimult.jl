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

    println(" Test picture invariance")
    ρ = rand(PauliSum{N})
    o = rand(PauliSum{N})
    ev_s = tr( o * (L*ρ))
    ev_h = tr((L'*o) * ρ)
    @test norm(ev_s - ev_h) < 1e-12

    println(" Test SubspaceDissipator")
    N = 4
    subspace = rand(PauliSum{N})
    D = SubspaceDissipator{N}(subspace, .1)
    display(D)
    display(subspace)
    display(collect(keys(subspace))[1])
    o = Pauli(collect(keys(subspace))[1])
    subspace += rand(PauliSum{N})
    display(o)
    display(subspace)
    println()
    display(D*subspace)
    return 
end

# test1()