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
    add_hamiltonian!(L, OpenSCI.heisenberg_1D(N, 1.1, 1.2, 1.3))
    add_channel_dephasing!(L, .1)
    add_channel_depolarizing!(L, .1)
    add_channel_amplitude_damping!(L, .1)
    display(L)



    eigvals, eigvecs = eigs(LinearMap(L), nev=5, which=:LR, tol=1e-5)
   
    for i in 1:length(eigvals)
        @printf(" %3i %12.8f %12.8fi\n", i, real(eigvals[i]), imag(eigvals[i]))
    end
    
    return L, eigvals, eigvecs
end

run()