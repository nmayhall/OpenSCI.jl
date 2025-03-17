using PauliOperators
using OpenSCI
using Test
using LinearAlgebra
using Printf
using Random
using Plots
using Statistics
using BenchmarkTools
# using DifferentialEquations
# using KrylovKit
# using Arpack

@inline commute(p1, p2) = iseven(count_ones(p1.x & p2.z) - count_ones(p1.z & p2.x)) 
function expectation_value(p::PauliBasis{N}, ket::Ket{N}) where N
    p.x == 0 || return 0.0
    
    count_ones(p.z & ket.v) % 2 == 0 || return -1
    return 1
end
expectation_value(p::Pauli, ket) = expectation_value(PauliBasis(p), ket) * coeff(p)

function PauliOperators.clip!(ps::PauliSum; thresh=1e-5)
    filter!(p->maximum(abs2.(p.second)) > thresh, ps)
    return ps 
end

function evolve_bfs(o::PauliSum{N}, g; thresh=-1) where N
    ot = deepcopy(o)
    for (i,gi) in enumerate(g)
        coeffi = coeff(gi)
        geni = PauliBasis(gi)

        sin_branch = PauliSum(N)
        for (oj,coeffj) in ot
            if commute(oj,gi) == false
                ot[oj] *= cos(coeffi)
                os = 1im * sin(coeffi) * oj * geni * coeffj
                sum!(sin_branch, os) 
                # sin_branch += 1im * sin(coeffi) * oj * geni * coeffj
            end
        end
        sum!(ot, clip!(sin_branch, thresh=thresh))
        PauliOperators.clip!(ot, thresh=thresh)
    end
    return ot
end

function form_U(o::PauliSum{N}, g) where N
    U = Matrix(1I, 2^N, 2^N)
    for (i,gi) in enumerate(g)
        U *= exp(-1im * Matrix(PauliBasis(gi)) * coeff(gi)/2)
    end
    return U
end
function evolve_exact(o::PauliSum{N}, g) where N
    U = form_U(o,g)
    return U' * Matrix(o) * U
end

function get_random_path(o, g)
    n_gen = length(g)
    sin_path = zeros(Bool, n_gen)
    cos_path = zeros(Bool, n_gen)

    ot = deepcopy(o)
    for (i,gi) in enumerate(g)
        if commute(ot,gi) == false
            # if tan(abs2(coeff(gi))) > rand()
            if rand() > .5
                sin_path[i] = 1
                ot = 1im*PauliBasis(gi)*ot
            else
                cos_path[i] = 1
            end
        end
    end
    return ot, sin_path, cos_path
end

function path_weight(sin_path::Vector{Bool}, cos_path::Vector{Bool}, g)
    length(sin_path) == length(cos_path) || throw(ErrorException)
    length(g) == length(cos_path) || throw(ErrorException)
    w = 1
    for i in 1:length(g)
        if cos_path[i] 
            w *= sin(coeff(g[i]))
        elseif sin_path[i] 
            w *= cos(coeff(g[i]))
        end
    end
    return w
end

function relative_weight(s1, c1, s2, c2, g)
    rw = 1
    for i in 1:length(g)
        coeff_i = coeff(g[i])
    
        rw *= tan(coeff_i)^(s2[i] & c1[i])
        rw *= cot(coeff_i)^(c2[i] & s1[i])
        rw *= cos(coeff_i)^(c2[i] & (s1[i] ⊽ c1[i]))
        rw *= sin(coeff_i)^(s2[i] & (s1[i] ⊽ c1[i]))
        rw *= sec(coeff_i)^(c1[i] & (s2[i] ⊽ c2[i]))
        rw *= csc(coeff_i)^(s1[i] & (s2[i] ⊽ c2[i]))

    end
    return rw
end

function monte_carlo(nsteps, o::Pauli{N}, g) where N
    trajectory = Vector{Pauli{N}}([]) 
    o_curr, s_curr, c_curr = get_random_path(o, g)
    push!(trajectory, o_curr)
    for step in 1:nsteps
        o_prop, s_prop, c_prop = get_random_path(o, g)
        rw = relative_weight(s_curr, c_curr, s_prop, c_prop, g)
        # w1 = path_weight(s_curr, c_curr, g)
        # w2 = path_weight(s_prop, c_prop, g)
        # rw = w2/w1
        # if abs(imag(rw)) > 1e-15
        #     throw(ErrorException)
        # end

        # if real(rw) < -1e-15
        #     @show rw
        #     throw(ErrorException)
        # end
        # tmp_curr = c_curr .- s_curr
        # tmp_prop = c_prop .- s_prop
        # display(o_curr)
        # display(o_prop)
        # @show abs(rw) tmp_curr tmp_prop 
        if abs2(rw)^2 > rand()
            # println("yes")
            o_curr = o_prop
            s_curr .= s_prop
            c_curr .= c_prop
        # else
            # println("no")
        end
        push!(trajectory, o_curr)
    end
    return trajectory
end

function run()
    N = 4
    α = 1
    k = 20 
    H = OpenSCI.heisenberg_1D(N, π/2, π/2, π/2, x=π/4*α)
    # display(H)

    Random.seed!(2)

    generators = Vector{Pauli{N}}([])
   
    for ki in 1:k
        for (i,j) in H
            push!(generators, i*j)
        end
    end
   
    o = Pauli(N,Z=[1])
    # for i in 2:N
    #     o += Pauli(N, Z=[i])
    # end
    ψ = Ket(N,0)

    # # display(O)
    # ot, sin_path, cos_path = get_path(o,generators)
    # display(ot) 
    # display(sin_path)
    # display(cos_path)

    # ot1, sin_path1, cos_path1 = get_path(o,generators)
    # ot2, sin_path2, cos_path2 = get_path(o,generators)
    # w1 = path_weight(sin_path1, cos_path1, generators)
    # w2 = path_weight(sin_path2, cos_path2, generators)
    # display(w1)
    # display(w2)
    # rw = relative_weight(sin_path1, cos_path1, sin_path2, cos_path2, generators)
    # # @printf(" Relative Weight direct = %12.8f %12.8fi\n", real(w1/w2), imag(w1/w2))
    # @printf(" Relative Weight        = %12.8f %12.8fi\n", real(rw), imag(rw))
    
    n_runs = 10
    data = zeros(Float64, n_runs)
    for i in 1:n_runs
        ev = get_stochastic_estimate(o, generators, ψ, n_samples=1000)
        # @printf(" Expectation Value mc    = %12.8f %12.8fi\n", real(ev), imag(ev))
        data[i] = ev
    end
   
    average = mean(data)
    @printf(" Expectation Value mc    = %12.8f ± %12.8f\n", mean(data), 3*std(data)/sqrt(n_runs))

    
    scatter(data, color=:grey, alpha = 0.1)
    plot!([s/i for (i,s) in enumerate(cumsum(data))], linewidth=2)
    savefig("./plot.pdf")

    # return
    # n_samples = 10000
    # ψ = Ket(N,0)
    # display(ψ)
    # traj = monte_carlo(n_samples, O, generators)
    # ev = 0
    # ot = PauliSum(N)
    # for i in traj
    #     ot += i
    #     evi = expectation_value(i,ψ) 
    #     # display(i)
    #     # display(evi)
    #     ev += evi/n_samples
    #     # @printf("%12.8f + %12.8fi  %20s\n", real(ev), imag(ev), i)
    # end
    # Z = 0
    # for (pi,ci) in ot
    #     Z += ci'*ci
    # end
    # @show Z
    # ot = ot * (1/sqrt(Z))
    # # ot *= 1/length(traj)
    # display(ot)
    # # ev = ev/length(traj)
    # # ev = ev/sqrt(N)
    # @printf(" Expectation Value mc    = %12.8f %12.8fi\n", real(ev), imag(ev))


    obfs = evolve_bfs(PauliSum(o), generators, thresh=1e-6)
    # display(obfs)
    @printf(" Number of BFS operators: %i\n", length(obfs))
    ev2 = 0
    for (p,coeff) in obfs
        ev2 += expectation_value(p,ψ) * coeff
    end
    @printf(" Expectation Value bfs   = %12.8f\n", ev2)
    
    omat = evolve_exact(PauliSum(o), generators)
    ev3 = Vector(ψ)' * omat * Vector(ψ) 
    @printf(" Expectation Value exact = %12.8f + %12.8fi\n", real(ev3), imag(ev3))
end


function get_stochastic_path(o, g, sin_vec)
    n_gen = length(g)
    path = zeros(Int8, n_gen)
    ot = deepcopy(o)
    for (i,gi) in enumerate(g)
        if commute(ot,gi) == false
            if sin_vec[i] > rand()
                path[i] = -1
                ot = 1im*PauliBasis(gi)*ot
            else
                path[i] = 1
            end
        end
    end
    return ot, path
end

function get_stochastic_estimate(o, g, ψ; n_samples=1000)
    # sin_vec = [sin(coeff(i))/(sin(coeff(i)) + cos(coeff(i))) for i in g]
    sin_vec = [sin(coeff(i))^2 for i in g]
    sin_vec = real.(sin_vec)
 
    n_gen = length(g)
    
    paths = Dict{NTuple{n_gen, Int8}, Float64}()

    for i in 1:n_samples
        o_i, path_i = get_stochastic_path(o, g, sin_vec)
        path_ii = ntuple(i->path_i[i], n_gen)
        if haskey(paths, path_ii)
            paths[path_ii] += 1 
        else
            paths[path_ii] = 1
        end
    end
   
    ev = 0
    for (path, Ni) in paths
        ot = deepcopy(o)
        for (i, branch) in enumerate(path)
            if branch == -1 
                ot = 1im*PauliBasis(g[i])*ot
            end
        end
        # prob = sqrt(Ni/n_samples)
        # ev += expectation_value(ot, ψ)*prob
        ev += expectation_value(ot, ψ)*Ni/n_samples
    end
    return ev
end

function test2()
    N = 2

    # Random.seed!(2)

    generators = Vector{Pauli{N}}([])
    α = π/4*.8
    k = 2 
    for i in 1:k
        push!(generators, Pauli("XX")*α)
        push!(generators, Pauli("YY")*α)
        push!(generators, Pauli("XI")*α)
        push!(generators, Pauli("YI")*α)
        push!(generators, Pauli("ZI")*α)
    end
    o = Pauli(N,X=[1])

    ψ = Ket(N,0)

    n_runs = 1000
    data = zeros(Float64, n_runs)
    for i in 1:n_runs
        ev = get_stochastic_estimate(o, generators, ψ, n_samples=1000)
        # @printf(" Expectation Value mc    = %12.8f %12.8fi\n", real(ev), imag(ev))
        data[i] = ev
    end
    
    average = mean(data)
    @printf(" Expectation Value mc    = %12.8f ± %12.8f\n", mean(data), 3*std(data)/sqrt(n_runs))

    
    scatter(data, color=:grey)
    plot!([s/i for (i,s) in enumerate(cumsum(data))])
    savefig("./plot.pdf")

    # Now exact
    obfs = evolve_bfs(PauliSum(o), generators)
    display(obfs)
    ev2 = 0
    for (p,coeff) in obfs
        ev2 += expectation_value(p,ψ) * coeff
    end
    @printf(" Expectation Value bfs   = %12.8f\n", ev2)
    
    omat = evolve_exact(PauliSum(o), generators)
    ev3 = Vector(ψ)' * omat * Vector(ψ) 
    @printf(" Expectation Value exact = %12.8f + %12.8fi\n", real(ev3), imag(ev3))
end

function test1()
    N = 5
    H = OpenSCI.heisenberg_1D(N, .2, 1.0, 1.2)
    # display(H)

    generators = Vector{Pauli{N}}([])
    k = 1
    for ki in 1:k
        for (i,j) in H
            push!(generators, i*j)
        end
    end
    

    O = PauliSum(Pauli(N,Z=[1]))

    O1 = evolve_exact(O, generators)
    O2 = evolve_bfs(O, generators)

    # println(" O(0)")
    # display(O)
    # println(" evolve: ")
    # println(" O(T)")
    # # display(O1)
    # display(O2)
    @test norm(O1 - Matrix(O2)) < 1e-12
end

run()
# test1()
# test2()