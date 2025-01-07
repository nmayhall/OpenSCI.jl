using PauliOperators
using OrderedCollections
# using LinearAlgebra
# using Printf
# using Dictionaries
# using StaticArrays

SparseDyadVectors{N,T} = OrderedDict{Dyad{N}, Vector{T}} 
function SparseDyadVectors(ds::DyadSum{N,T}) where {N,T}
    sdv = OrderedDict{Dyad{N}, Vector{T}}()
    for (dyad,coeff) in ds
        sdv[dyad] = [coeff] 
    end
    return sdv
end

function clip!(ps::SparseDyadVectors{N,T}; thresh=1e-16) where {N,T}
    filter!(p->maximum(abs.(p.second)) > thresh, ps.ops)
end

"""
    H::PauliSum{N}
    L::Vector{PauliSum{N}}
    γ::Vector{Float64}

This is general enough to realize infinitesmal evolution of any CPT map. 
`H`: generates unitary evolution 
`L`: Is a vector of Jump operators
`γ`: Is the vector of Jump operator coefficients (diagonal form assumed)
"""
struct Lindbladian{N}
    H::PauliSum{N}          
    L::Vector{PauliSum{N}}
    γ::Vector{Float64}
end

function Lindbladian(N)
    return Lindbladian{N}(PauliSum(N), Vector{PauliSum{N}}([]), Vector{Float64}([]))
end

function Lindbladian(H, L, γ)
    return Lindbladian{N}(H,L,γ)
end

function Base.rand(T::Type{Lindbladian{N}}; nH=2, nL=2) where N
    Λ = Lindbladian(N)
    for i in 1:nH
        opi = rand(ScaledPauli{N})
        opi = opi + opi'
        # Λ.H += opi
        sum!(Λ.H,opi)
    end
    for i in 1:nL
        opi = rand(Pauli{N})
        # for j in 1:nL
        #     opj = rand(ScaledPauli{N})
        #     # is_diagonal(opj) == false || continue
        #     opi += opj  
        # end
        # opi = opi + opi'
        push!(Λ.L,PauliSum(opi))
        push!(Λ.γ,rand())
    end
    return Λ
end

function Base.string(Λ::Lindbladian{N}) where N
    out = "H:\n"
    for (p,c) in Λ.H.ops
        out *= @sprintf(" % 12.8f + % 12.8fi %s\n", real(c), imag(c), p)
    end
    out *= "\n"
    out *= "L:\n"
    for i in 1:length(Λ.L)
        out *= @sprintf("%4i : γ = % 12.8f + % 12.8fi\n", i, real(Λ.γ[i]), imag(Λ.γ[i]))
        for (p,c) in Λ.L[i].ops
            out *= @sprintf("   % 12.8f + % 12.8fi %s\n", real(c), imag(c), p)
        end
        # out *= string(Λ.γ[i]) * " : " * string(Λ.L[i]) * "\n"
    end
    return out 
end

Base.display(Λ::Lindbladian) = println(string(Λ))
function Base.display(sdv::SparseDyadVectors{N,T})  where {N,T}
    idx = 1
    for (d,coeffs) in sdv
        @printf(" %6i %12s", idx, d)
        for coeff in coeffs
            @printf(" % 10.8f % 10.8fi", real(coeff), imag(coeff))
        end
        @printf("\n")
        idx += 1
    end
end

Base.Vector(sdv::SparseDyadVectors{N,T}; state=1)  where {N,T} = [i[state] for i in values(sdv)]
function Base.Matrix(sdv::SparseDyadVectors{N,T})  where {N,T}
    ni, nj = size(sdv)
    out = zeros(T, ni, nj)
    i = 1
    for (d,coeffs) in sdv
        for j in 1:nj 
            out[i,j] = coeffs[j]
        end
        i += 1
    end
    return out
end

function Base.size(sdv::SparseDyadVectors)
    ni = length(sdv)
    for (d,c) in sdv
        return (ni, length(c)) 
    end
end

function Base.sum!(sdv::SparseDyadVectors{N,T}, dyad::Dyad{N}, coeffs::Vector{T}) where {N,T}
    if haskey(sdv, dyad)
        sdv[dyad] .+= coeffs
    else
        sdv[dyad] = coeffs
    end
end

function todense(sdv::SparseDyadVectors{N,T}) where {N,T}
    out = zeros(T, 4^N, size(sdv)[2])
    for (d,c) in sdv
        out[vec(d),:] .= c
    end
    return out
end