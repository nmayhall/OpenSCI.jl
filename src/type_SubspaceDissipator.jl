"""
    set::PauliSum{N}
    γ::Float64

Dampen operators outside of operators in `set`. 
As this is a diagonal operator, exponentiation is trivial.
D = -γ∑_{i∉set} |Qᵢ)(Qᵢ|
D = -γQ = -γ(I-P)
"""
struct SubspaceDissipator{N} <: SuperOperator{N}
    set::PauliSum{N, ComplexF64}
    γ::Float64
end

function Base.:*(D::SubspaceDissipator, o::PauliSum)
    out = deepcopy(o)
    return lmul!(D,out)
end

Base.:*(D::SubspaceDissipator, o::Union{Pauli,PauliBasis}) = D*PauliSum(o)

"""
    LinearAlgebra.lmul!(D::SubspaceDissipator, o::PauliSum)

Peform D*o and overwrite o
"""
function LinearAlgebra.lmul!(D::SubspaceDissipator, o::PauliSum)
    for (oi,coeff) in o
        if haskey(D.set, oi) == false
            o[oi] *= -D.γ
        end
    end
    return o
end
