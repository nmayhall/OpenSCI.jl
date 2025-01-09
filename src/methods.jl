using PauliOperators
using OpenSCI


function Base.:*(L::Lindbladian{N}, ρ::DyadSum{N,T}) where {N,T}
    # out = DyadSum(N,T=T)
    # rmul!(out, L, ρ)
    dρ = -1im * (L.H*ρ - ρ*L.H)
    for i in 1:length(L.γ)
        Li = L.L[i]
        dρ += L.γ[i] * (Li * ρ * Li')
        LL = Li' * Li
        dρ -= 0.5*L.γ[i]*(LL * ρ + ρ * LL)
    end
    return dρ 
end
function Base.:*(L::Lindbladian{N}, ρd::Dyad{N}, T=ComplexF64) where N
    return L*DyadSum(ρd, T=T) 
end

Base.vec(d::Dyad{N}) where N = 1 + d.ket.v + d.bra.v*(BigInt(2)^N)

function Base.Matrix(L::Lindbladian{N}) where N
    Lmat = zeros(ComplexF64, 4^N, 4^N)
    for kr in 0:2^N-1
        for br in 0:2^N-1
            dyad_r = Dyad(N, kr, br)
            
            σ = L * dyad_r        
            for (dyad_l, coeff) in σ
                Lmat[vec(dyad_l), vec(dyad_r)] = coeff
            end
        end
    end
    return Lmat
end

function add_hamiltonian!(L::Lindbladian{N}, H::PauliSum{N}) where N
    sum!(L.H, H)
    return L
end

function add_channel_amplitude_damping!(L::Lindbladian{N}, γ::Real) where N
    return add_channel_amplitude_damping!(L, [γ for i in 1:N])
end

function add_channel_dephasing!(L::Lindbladian{N}, γ::Real) where N
    return add_channel_dephasing!(L, [γ for i in 1:N])
end

function add_channel_depolarizing!(L::Lindbladian{N}, γ::Real) where N
    return add_channel_depolarizing!(L, [γ for i in 1:N])
end


function add_channel_amplitude_damping!(L::Lindbladian{N}, γ::Vector{<:Real}) where N
    length(γ) == N || throw(ArgumentError("Length of γ must be equal to N"))
    for i in 1:N
        pi = PauliSum(N)
        pi += 1/sqrt(2)*(Pauli(N, X=[i]) + 1im * Pauli(N, Y=[i]))
        push!(L.L, pi)
        push!(L.γ, γ[i])
    end
end


function add_channel_dephasing!(L::Lindbladian{N}, γ::Vector{<:Real}) where N
    length(γ) == N || throw(ArgumentError("Length of γ must be equal to N"))
    for i in 1:N
        pi = PauliSum(N)
        pi += 1*Pauli(N, Z=[i])
        push!(L.L, pi)
        push!(L.γ, γ[i])
    end
end


function add_channel_depolarizing!(L::Lindbladian{N}, γ::Vector{<:Real}) where N
    length(γ) == N || throw(ArgumentError("Length of γ must be equal to N"))
    for i in 1:N
        pi = PauliSum(N)
        pi += Pauli(N, X=[i])
        push!(L.L, pi)
        push!(L.γ, γ[i])
        
        pi = PauliSum(N)
        pi += Pauli(N, Y=[i])
        push!(L.L, pi)
        push!(L.γ, γ[i])
        
        pi = PauliSum(N)
        pi += Pauli(N, Z=[i])
        push!(L.L, pi)
        push!(L.γ, γ[i])
    end
end

