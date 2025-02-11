"""
    H::PauliSum{N}
    L::Vector{PauliSum{N}}
    γ::Vector{Float64}

This is general enough to realize infinitesmal evolution of any CPT map. 
`H`: generates unitary evolution 
`L`: Is a vector of Jump operators
`γ`: Is the vector of Jump operator coefficients (diagonal form assumed)
"""
struct Lindbladian{N} <: SuperOperator{N}
    H::PauliSum{N}          
    L::Vector{PauliSum{N}}
    γ::Vector{Float64}
end

function Lindbladian(N)
    return Lindbladian{N}(PauliSum(N, ComplexF64), Vector{PauliSum{N}}([]), Vector{Float64}([]))
end

function Lindbladian(H, L, γ)
    return Lindbladian{N}(H,L,γ)
end

Base.adjoint(d::Lindbladian{N}) where N = Adjoint(d)
Base.parent(d::Adjoint{<:Any, <:Lindbladian}) = d.parent

function Base.rand(T::Type{Lindbladian{N}}; nH=2, nL=2) where N
    Λ = Lindbladian(N)
    for i in 1:nH
        opi = rand(Pauli{N})
        opi = opi + opi'
        # Λ.H += opi
        sum!(Λ.H,opi*rand())
    end
    for i in 1:nL
        opi = rand(Pauli{N})
        push!(Λ.L,PauliSum(opi))
        push!(Λ.γ,rand())
    end
    return Λ
end

function Base.string(Λ::Lindbladian{N}) where N
    out = "H:\n"
    for (p,c) in Λ.H
        out *= @sprintf(" % 12.8f + % 12.8fi %s\n", real(c), imag(c), p)
    end
    out *= "\n"
    out *= "L:\n"
    for i in 1:length(Λ.L)
        out *= @sprintf("%4i : γ = % 12.8f + % 12.8fi\n", i, real(Λ.γ[i]), imag(Λ.γ[i]))
        for (p,c) in Λ.L[i]
            out *= @sprintf("   % 12.8f + % 12.8fi %s\n", real(c), imag(c), p)
        end
        # out *= string(Λ.γ[i]) * " : " * string(Λ.L[i]) * "\n"
    end
    return out 
end

Base.display(Λ::Lindbladian) = println(string(Λ))


"""
    LinearAlgebra.mul!(σ::Vector{T}, L::Lindbladian{N}, v::Vector{T}) where {N,T}

Compute the action of `L*v` and store in σ. 
# Arguments
- `σ::Vector{T}`: output vector (given in Dyad Basis)
- `L::Lindbladian{N}`: Lindbladian SuperOperator
- `v::Vector{T}`: input vector (assumed to be coefficients in the Dyad Basis) 
"""
function LinearAlgebra.mul!(σ::Vector{T}, L::Lindbladian{N}, v::Vector{T}) where {N,T}
    for rdyadbasis in DyadBasis{N}
        rcoeff = v[index(rdyadbasis)]

        # Unitary part
        for (pauli,coeff_h) in L.H
            ldyad = -1im * coeff_h * (pauli * rdyadbasis)
            lcoeff = coeff(ldyad)
            ldyadbasis = DyadBasis(ldyad)
            σ[index(ldyadbasis)] += lcoeff
            
            ldyad = -1im * coeff_h * (rdyadbasis * pauli)
            lcoeff = coeff(ldyad)
            ldyadbasis = DyadBasis(ldyad)
            σ[index(ldyadbasis)] -= lcoeff
        end

        # Non-unitary part
        for i in 1:length(L.γ)
            coeff_i = L.γ[i] 
            for (pauli_j, coeff_j) in L.L[i]
                for (pauli_k, coeff_k) in L.L[i]
                    # Li v Li'
                    ldyad = rcoeff * coeff_i * coeff_j * coeff_k' * (pauli_j * (rdyadbasis * pauli_k'))
                    lcoeff = coeff(ldyad)
                    ldyadbasis = DyadBasis(ldyad)
                    σ[index(ldyadbasis)] += lcoeff
                    
                    # -1/2 Li'Li v
                    ldyad = (-1/2) * rcoeff * coeff_i * coeff_j' * coeff_k * (pauli_j' * pauli_k * rdyadbasis)
                    lcoeff = coeff(ldyad)
                    ldyadbasis = DyadBasis(ldyad)
                    σ[index(ldyadbasis)] += lcoeff
                    
                    # -1/2 v Li'Li
                    ldyad = (-1/2) * rcoeff * coeff_i * coeff_j' * coeff_k * (rdyadbasis * pauli_j' * pauli_k)
                    lcoeff = coeff(ldyad)
                    ldyadbasis = DyadBasis(ldyad)
                    σ[index(ldyadbasis)] += lcoeff
                end
            end
        end
    end

    # # Unitary part
    # for (rdyad, rcoeffs) in ρ

    #     σi = -1im * (L.H * rdyad - rdyad * L.H)
        
    #     for i in 1:length(L.γ)
    #         Li = L.L[i]
    #         LL = Li' * Li
           
    #         σi +=  L.γ[i] * (Li * rdyad * Li')
    #         σi -= 0.5*L.γ[i]*(LL*rdyad + rdyad*LL)
    #     end
        
    #     for (ldyad,lcoeff) in σi 
    #         sum!(σ, ldyad, lcoeff .* rcoeffs )
    #     end
    # end
    # return σ 
end
