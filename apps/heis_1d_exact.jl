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
using LinearMaps

# @testset "tests1.jl" begin
function run()
    N = 6
    L = Lindbladian(N)
    add_hamiltonian!(L, OpenSCI.heisenberg_1D(N, 1.0, 1.0, 1.0))
    add_channel_dephasing!(L, .1)
    add_channel_depolarizing!(L, .0)
    add_channel_amplitude_damping!(L, .0)
    display(L)

    R = 3

    Lmap = LinearMap(L)
    time = @elapsed e0, v, info = KrylovKit.eigsolve(Lmap, R)

    # time = @elapsed e0, v, info = KrylovKit.eigsolve(Lmap, R, :LR,
    #     verbosity=1,
    #     maxiter=200,
    #     #krylovdim=20, 
    #     issymmetric=false,
    #     ishermitian=false,
    #     tol=1e-6)


    # eigvals, eigvecs = eigsolve(Lmap, nev=5, which=:LR, tol=1e-5)
    # eigvals, eigvecs = eigs(Lmap, nev=5, which=:LR, tol=1e-5)
   
    for i in 1:length(eigvals)
        @printf(" %3i %12.8f %12.8fi\n", i, real(eigvals[i]), imag(eigvals[i]))
    end
    
    return L, eigvals, eigvecs
end

run()