abstract type AbstractKet{N} end


struct PauliKet{N} <: AbstractKet{N}
    operator::FixedPhasePauli{N}
end
struct DyadKet{N} <: AbstractKet{N}
    operator::Dyad{N}
end

Base.adjoint(d::PauliKet{N}) where N = Adjoint(d)
Base.adjoint(d::DyadKet{N}) where N = Adjoint(d)
Base.parent(d::Adjoint{<:Any, <:PauliKet}) = d.parent
Base.parent(d::Adjoint{<:Any, <:DyadKet}) = d.parent


"""
    Base.adjoint(p::FixedPhasePauli{N}) where N

Here we assume that the pauli's in the Ket are Hermitian
"""
function Base.adjoint(p::FixedPhasePauli{N}) where N
    return p
end

PauliOperators.FixedPhasePauli(p::FixedPhasePauli) = p 
PauliOperators.FixedPhasePauli(p::Pauli) = p.pauli 
PauliOperators.FixedPhasePauli(p::ScaledPauli) = p.pauli 
PauliOperators.Dyad(d::Dyad) = d
PauliOperators.Dyad(d::ScaledDyad) = d.dyad

Ket(p::Union{FixedPhasePauli{N}, Pauli{N}, ScaledPauli{N}}) where N = PauliKet{N}(FixedPhasePauli(p))
Ket(d::Union{Dyad{N}, ScaledDyad{N,T}}) where {N,T} = DyadKet{N}(Dyad(d))

Base.:*(a::Adjoint{<:Any, <:AbstractKet{N}}, b::AbstractKet{N}) where N = dot(a.parent.operator',b.operator)

function LinearAlgebra.dot(a::FixedPhasePauli{N}, b::Dyad{N}) where N
    i1 = iseven(count_ones(a.z & b.bra.v)) 
    i2 = b.bra.v == (a.x ⊻ b.ket.v)
    return -(-1)^i1*i2 
end
function LinearAlgebra.dot(b::Dyad{N}, a::FixedPhasePauli{N}) where N
    i1 = iseven(count_ones(a.z & b.bra.v)) 
    i2 = b.bra.v == (a.x ⊻ b.ket.v)
    return -(-1)^i1*i2 
end
