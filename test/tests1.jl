using PauliOperators
using OpenSCI
using Test
using LinearAlgebra
using Printf
using Random

# @testset "tests1.jl" begin
function run1()
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
end

function run2()

    #### Now with non-unitary part 
    N = 2
    Λ = rand(Lindbladian{N}, nH=100, nL=1)
    display(Λ)
    @show ρ = DyadSum(Dyad(N,0,0)) 

    err = Matrix(Λ)*vec(Matrix(ρ)) - vec(Matrix(Λ*ρ))
    @test isapprox(norm(err), 0, atol=1e-14)

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
    ρ = Matrix{Float64}(I, 2^N, 2^N)
    
    F = svd(rand(2^N, 2^N))
    ρ = F.U * Diagonal(F.S) * F.U' / sum(F.S)

    println(" Eigenvalues of ρ ", eigvals(ρ))
    println("\n Trace of ρ: ", tr(ρ))
    println(" ρ")
    display(ρ)
    Λρ = Matrix(Λ)*vec(ρ)
    Λρ = reshape(Λρ, (2^N, 2^N))

    println(" Λρ")
    display(Λρ)
    println("\n Trace of Λρ: ", tr(Λρ))
    println(" Eigenvalues of ρ ", eigvals(Λρ))

    println("\n  Do my jump operators sum to I?")
    tmp = PauliSum(N) 
    for Li in Λ.L
        tmp += Li*Li'
    end
    display(tmp)
end

run1()
run2()
