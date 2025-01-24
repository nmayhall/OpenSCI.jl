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

LinearAlgebra.dot(a::FixedPhasePauli{N}, b::FixedPhasePauli{N}) where N = a == b
LinearAlgebra.dot(a::Dyad{N}, b::Dyad{N}) where N = (a.ket.v == b.bra.v) && (a.bra.v == b.ket.v)


PauliTypes{N} = Union{Pauli{N}, ScaledPauli{N}}

PauliOperators.get_coeff(p::FixedPhasePauli) = 1
PauliOperators.get_coeff(p::Dyad) = 1
PauliOperators.get_coeff(p::ScaledDyad) = p.coeff
PauliOperators.FixedPhasePauli(p::FixedPhasePauli) = p 
PauliOperators.FixedPhasePauli(p::Pauli) = p.pauli 
PauliOperators.FixedPhasePauli(p::ScaledPauli) = p.pauli 
PauliOperators.Dyad(d::Dyad) = d
PauliOperators.Dyad(d::ScaledDyad) = d.dyad

LinearAlgebra.dot(a::PauliTypes, b::PauliTypes) = get_coeff(a)*get_coeff(b) * dot(FixedPhasePauli(a), FixedPhasePauli(b)) 
LinearAlgebra.dot(a::PauliTypes, b::FixedPhasePauli) = get_coeff(a)*get_coeff(b) * dot(FixedPhasePauli(a), FixedPhasePauli(b)) 
LinearAlgebra.dot(a::FixedPhasePauli, b::PauliTypes) = get_coeff(a)*get_coeff(b) * dot(FixedPhasePauli(a), FixedPhasePauli(b)) 


LinearAlgebra.dot(a::ScaledDyad, b::ScaledDyad) = get_coeff(a)*get_coeff(b) * dot(Dyad(a), Dyad(b)) 
LinearAlgebra.dot(a::ScaledDyad, b::Dyad) = get_coeff(a)*get_coeff(b) * dot(Dyad(a), Dyad(b)) 
LinearAlgebra.dot(a::Dyad, b::ScaledDyad) = get_coeff(a)*get_coeff(b) * dot(Dyad(a), Dyad(b)) 

LinearAlgebra.dot(a::FixedPhasePauli, b::Dyad) = get_coeff(a)*get_coeff(b) * dot(FixedPhasePauli(a), Dyad(b)) 
LinearAlgebra.dot(a::FixedPhasePauli, b::ScaledDyad) = get_coeff(a)*get_coeff(b) * dot(FixedPhasePauli(a), Dyad(b)) 
LinearAlgebra.dot(a::Dyad, b::FixedPhasePauli) = get_coeff(a)*get_coeff(b) * dot(Dyad(a), FixedPhasePauli(b)) 
LinearAlgebra.dot(a::ScaledDyad, b::FixedPhasePauli) = get_coeff(a)*get_coeff(b) * dot(Dyad(a), FixedPhasePauli(b)) 

LinearAlgebra.dot(a::PauliTypes, b::Dyad) = get_coeff(a)*get_coeff(b) * dot(FixedPhasePauli(a), Dyad(b)) 
LinearAlgebra.dot(a::PauliTypes{N}, b::ScaledDyad{N,T}) where {N,T} = get_coeff(a)*get_coeff(b) * dot(FixedPhasePauli(a), Dyad(b)) 
LinearAlgebra.dot(a::Dyad, b::PauliTypes) = get_coeff(a)*get_coeff(b) * dot(Dyad(a), FixedPhasePauli(b)) 
LinearAlgebra.dot(a::ScaledDyad, b::PauliTypes) = get_coeff(a)*get_coeff(b) * dot(Dyad(a), FixedPhasePauli(b)) 

