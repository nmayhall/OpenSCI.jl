using PauliOperators
using OpenSCI
using Test
using LinearAlgebra
using Printf
using Random
using Plots
using DifferentialEquations

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
    
    err = -1im*(vec(Matrix(Λ.H)*Matrix(ρ) - Matrix(ρ)*Matrix(Λ.H))) - vec(Matrix(Λ*ρ))
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
    N = 1
    Λ = rand(Lindbladian{N}, nH=100, nL=0)
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
    println("\n Trace of ρ+Λρ: ", tr(ρ + Λρ))
    println(" Eigenvalues of ρ+Λρ ", eigvals(ρ + Λρ))

    println("\n  Do my jump operators sum to I?")
    tmp = PauliSum(N) 
    for Li in Λ.L
        tmp += Li*Li'
    end
    display(tmp)


    println("\n Let's try again")
    display(ρ)
    display(ρ + -1im*(Matrix(Λ.H)*ρ - ρ*Matrix(Λ.H)))
end

function run3()
    N = 3
    dim = 2^N
    # Λ = rand(Lindbladian{N}, nH=100, nL=10)
    Λ = Lindbladian(N)
    add_hamiltonian!(Λ, OpenSCI.heisenberg_1D(N, 10.1, 20.2, 2.3))
    add_channel_dephasing!(Λ, .1)
    add_channel_depolarizing!(Λ, .1)
    # add_channel_amplitude_damping!(Λ, 1.)
    display(Λ)
    
    println("\n EigenValues of Λ")
    L_evals = eigvals(Matrix(Λ))
    for i in 1:length(L_evals) 
        @printf(" %4i %12.8f %12.8fi\n", i, real(L_evals[i]), imag(L_evals[i]))
    end
    
    U,s,V = svd(Matrix(Λ))
    println("\n Singular Values of Λ")
    for i in 1:length(s)
        @printf(" %4i %12.8f %12.8fi\n", i, real(s[i]), imag(s[i]))
    end
    V = V ./ sqrt(2^N)
    println("\n ρss: ")
    ρss = reshape(V[:,end], (2^N, 2^N))
    display(real(ρss))

    println("\n ρss-1: ")
    display(real(reshape(V[:,end-1], (2^N, 2^N))))

    @show ρ = DyadSum(Dyad(N,0,0))


    Lmat = Matrix(Λ)*(1+0.0im)
    ρmat = vec(Matrix(ρ))*(1+0.0im)
   
    F = eigen(Lmat)
    println("\n Eigenvalues of Λ")
    for i in 1:length(s)
        @printf(" %4i %12.8f %12.8fi\n", i, real(F.values[i]), imag(F.values[i]))
    end
    println(" Trace of ρss: ", tr(reshape(F.vectors[:,end], (2^N, 2^N)))/sqrt(2^N))
    println("")


    f(du, u, p, t) = du .= Lmat*u
    tspan = (0.0, 1.0)
    prob = ODEProblem(f, ρmat, tspan) 
    # f(u, p, t) = Λ*u 
    # tspan = (0.0, 1.0)
    # prob = ODEProblem(f, ρ, tspan) 
    
    sol = solve(prob, reltol = 1e-8, abstol = 1e-8, saveat = 0.01)

    gs_prob = [abs(i[1]*i[1]') for i in sol.u]
    t = [i for i in sol.t]
    plot(t, gs_prob, label = "gs_prob")
    for ei in L_evals[1:5] 
        println(" ei = ", ei)
        plot!(t, [abs(exp(ei*ti)) for ti in t], linestyle = :dash, label = "exp(ei*ti)", legend=false)
    end

    for i in 1:length(sol.t)
        ρt = reshape(sol.u[i], (2^N, 2^N))
        @printf(" tr(ρt) = %12.8f %12.8fi\n", real(tr(ρt)), imag(tr(ρt)))
    end
    savefig("./plot.pdf")
    
    ρT = reshape(sol.u[end], (2^N, 2^N))
    println("\n ρT")
    display(ρT)
    println("\n eigvals(ρT)")
    display(eigvals(ρT))
    println("\n tr(ρT)")
    display(tr(ρT))

    println("\n tr(ρss*ρT)")
    display(tr(ρss*ρT))
end

run1()
run2()
run3()
