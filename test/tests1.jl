using PauliOperators
using OpenSCI
using Test
using LinearAlgebra

@testset "tests1.jl" begin

    N = 6

    @test vec(Dyad(N, 4, 8)) == 517 
    Λ = rand(Lindbladian{N})
    display(Λ)

    ρ = DyadSum(Dyad(N,0,0)) 
    display(ρ)

    σ = DyadSum(N)
    
    display(typeof(ρ))
    display(typeof(σ))
    # rmul!(σ, Λ, ρ)
    σ = Λ*ρ

    # test this with Matrix
    ρv = vec(Matrix(ρ))
    Lv = Matrix(Λ)
    
    err = Matrix(Λ)*vec(Matrix(ρ)) - vec(Matrix(Λ*ρ))
    @test isapprox(norm(err), 0, atol=1e-14)
    
end
