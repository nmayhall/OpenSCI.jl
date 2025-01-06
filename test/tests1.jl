using PauliOperators
using OpenSCI
using Test
using LinearAlgebra
using Printf
using Random

# @testset "tests1.jl" begin
function run()
    N = 2

    @test vec(Dyad(6, 4, 8)) == 517 
    Λ = rand(Lindbladian{N}, nH=100, nL=0)
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
    
    err = vec(Matrix(Λ.H)*Matrix(ρ) - Matrix(ρ)*Matrix(Λ.H)) - vec(Matrix(Λ*ρ))
    @test isapprox(norm(err), 0, atol=1e-14)

    U,s,V = svd(Lv)
    for i in 1:length(s)
        @printf(" %4i %12.8f %12.8fi\n", i, real(s[i]), imag(s[i]))
    end

   
    

    Λ = rand(Lindbladian{N}, nH=100, nL=0)

    display(Λ)
    F = eigen(Matrix(Λ.H))
    ρ = zeros(2^N, 2^N)
    for i in 1:length(F.values)
        ρ += rand()*F.vectors * F.vectors'
    end
    ρ = ρ/tr(ρ)
    # display(ρ)

    @show norm(Matrix(Λ)*vec(Matrix(ρ)))
   
    # return 

    #### Now with non-unitary part 
    N = 2
    Λ = rand(Lindbladian{N}, nH=100, nL=10)
    display(Λ)
    @show ρ = DyadSum(Dyad(N,0,0)) 

    err = Matrix(Λ)*vec(Matrix(ρ)) - vec(Matrix(Λ*ρ))
    @test isapprox(norm(err), 0, atol=1e-14)

    println("\n Singular Values of Λ")
    U,s,V = svd(Matrix(Λ))
    
    println("\n Singular Values of Λ")
    for i in 1:length(s)
        @printf(" %4i %12.8f %12.8fi\n", i, real(s[i]), imag(s[i]))
    end
    # display(V[:,4])
    ρss = real(reshape(V[:,4], (2^N, 2^N)))
    ρss = ρss/tr(ρss)
    
    println("\n ρss")
    display(ρss)
    @printf(" tr(ρss) = %12.8f %12.8fi\n", real(tr(ρss)), imag(tr(ρss)))

    Random.seed!(2)
    tmp = Matrix{Float64}(I, 2^N, 2^N)
    
    F = svd(rand(2^N, 2^N))
    tmp = F.U * Diagonal(F.S) * F.U'

    tmp /= tr(tmp) 
    println(" positive? ", eigvals(tmp))
    println("\n Trace of ρ")
    display(tr(tmp))
    display(tmp)
    tmp = Matrix(Λ)*vec(tmp)
    tmp = reshape(tmp, (2^N, 2^N))
    println("\n Trace of Λρ")
    display(tr(tmp))
    display(tmp)

    println("\n  Do my jump operators sum to I?")
    tmp = PauliSum(N) 
    for Li in Λ.L
        tmp += Li*Li'
    end
    display(tmp)
end
run()
