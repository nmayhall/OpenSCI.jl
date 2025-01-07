using PauliOperators
using OpenSCI
using Test
using LinearAlgebra
using Printf
using Random
using Plots
using DifferentialEquations

# @testset "tests1.jl" begin
function test1()
    N = 4
    L = rand(Lindbladian{N}, nH=10, nL=5)
    
    v0 = SparseDyadVectors(DyadSum(Dyad(N,0,0)))

    display(L)
    display(v0)
    display(Vector(v0))
    display(Matrix(v0))

    tmp1 = SparseDyadVectors{N,ComplexF64}()
    for i in 1:20
        tmp1[rand(Dyad{N})] = [randn() + 1im*randn() for j in 1:4]
    end
    @test norm(Matrix(tmp1) - Matrix(fill!(tmp1,Matrix(tmp1)))) < 1e-14

    println()
    display(L*v0)

    err = Matrix(L)*todense(v0) - todense(L*v0)
    @test isapprox(norm(err), 0, atol=1e-14)

    selected_ci(L, v0, max_iter_outer=5)

    # σ = multiply(L, v0)
    # Lmat = build_subspace_L(L, σ)
  
    Lmat = Matrix(L)
    _,s,V = svd(Lmat)

    @printf("\n Singular values of Lmat:\n")
    for i in 1:length(s)
        @printf(" %4i % 12.8f\n", i, s[i])
    end
end

test1()