using PauliOperators
using OpenSCI
using Test
using LinearAlgebra
using Printf
using Random
using Plots
using DifferentialEquations
using KrylovKit
using Arpack

# @testset "tests1.jl" begin
function run()
    N = 6
    L = Lindbladian(N)
    add_hamiltonian!(L, OpenSCI.heisenberg_1D(N, 1.0, 1.0, 1.0))
    add_channel_dephasing!(L, .1)
    add_channel_depolarizing!(L, .1)
    add_channel_amplitude_damping!(L, .1)
    display(L)



    # v0 = SparseDyadVectors(DyadSum(Dyad(N,0,0)), N=2)
    v0 = DyadBasis(N,0,0)
    v0 += DyadBasis(N, 0,1)
    v0 += DyadBasis(N, 1,0)
    v0 += DyadBasis(N, 1,1)
    v0 += DyadBasis(N, 2,2)
    v0 += DyadBasis(N, 2,0)
    v0 += DyadBasis(N, 0,2)
    
    v0 = SparseDyadVectors(v0, R=4)
    OpenSCI.eye!(v0)
    display(v0)
    
    v1 = selected_ci(L, v0, max_iter_outer=9, ϵsearch=1e-4, ϵdiscard=1e-6, thresh_conv=1e-4)

    # Lmat = Matrix(L)
    # l = eigvals(Lmat)
    # @printf("\n Exact Eigenvalues of Lmat:\n")
    # for i in length(l)-5:length(l)
    #     @printf(" %4i % 12.8f % 12.8fi\n", i, real(l[i]), imag(l[i]))
    # end

    return v1
end

run()