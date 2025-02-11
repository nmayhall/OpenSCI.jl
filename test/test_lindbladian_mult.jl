using PauliOperators
using OpenSCI
using Test
using LinearAlgebra
using Printf
using Random
using Plots

@testset "Lindbladian mul!" begin
    Random.seed!(1234)
    N = 3

    L = Lindbladian(N)
    add_hamiltonian!(L, OpenSCI.heisenberg_1D(N, 1.1, 1.2, 1.3))
    add_channel_dephasing!(L, .1)
    add_channel_depolarizing!(L, .1)
    add_channel_amplitude_damping!(L, .1)
    display(L)

    v0 = rand(DyadSum{N}, n_terms=100)
    v0 = v0 * (1/tr(v0))

    v = zeros(ComplexF64, 4^N)
    for (d,c) in v0
        v[index(d)] = c
    end
    v = v / norm(v)

    Lmat = Matrix(L)
    σref = Lmat * v 

    σtst = zeros(ComplexF64,size(σref)...)
    
    mul!(σtst, L, v)

    @test norm(σtst - σref) < 1e-15
end