function LinearAlgebra.dot(a::FixedPhasePauli{N}, b::Dyad{N}) where N
    i1 = iseven(count_ones(a.z & b.bra.v)) 
    i2 = b.bra.v == (a.x ⊻ b.ket.v)
    return -(-1)^i1*i2 
end

LinearAlgebra.dot(a::FixedPhasePauli{N}, b::FixedPhasePauli{N}) where N = a == b
LinearAlgebra.dot(a::Dyad{N}, b::Dyad{N}) where N = (a.ket == b.bra) && (a.bra == b.ket)

PauliTypes{N} = Union{Pauli{N}, ScaledPauli{N}}
DyadTypes{N} = Union{Dyad{N}, ScaledDyad{N}}

PauliOperators.get_coeff(p::FixedPhasePauli) = 1
PauliOperators.get_coeff(p::Dyad) = 1
PauliOperators.get_coeff(p::ScaledDyad) = p.coeff
# PauliOperators.FixedPhasePauli(p::FixedPhasePauli) = p 
PauliOperators.FixedPhasePauli(p::Pauli) = p.pauli 
PauliOperators.FixedPhasePauli(p::ScaledPauli) = p.pauli 
# PauliOperators.Dyad(d::Dyad) = d
PauliOperators.Dyad(d::ScaledDyad) = d.dyad

LinearAlgebra.dot(a::PauliTypes{N}, b::DyadTypes{N}) where N = get_coeff(a)*get_coeff(b) * dot(FixedPhasePauli(a), Dyad(b)) 
LinearAlgebra.dot(a::DyadTypes{N}, b::DyadTypes{N}) where N = get_coeff(a)*get_coeff(b) * dot(Dyad(a), Dyad(b)) 
LinearAlgebra.dot(a::DyadTypes{N}, b::PauliTypes{N}) where N = get_coeff(a)*get_coeff(b) * dot(Dyad(a), FixedPhasePauli(b)) 
LinearAlgebra.dot(a::PauliTypes{N}, b::PauliTypes{N}) where N = get_coeff(a)*get_coeff(b) * dot(FixedPhasePauli(a), FixedPhasePauli(b)) 
LinearAlgebra.dot(a::FixedPhasePauli{N}, b::PauliTypes{N}) where N = get_coeff(a)*get_coeff(b) * dot(a, FixedPhasePauli(b)) 
LinearAlgebra.dot(a::FixedPhasePauli{N}, b::DyadTypes{N}) where N = get_coeff(a)*get_coeff(b) * dot(a, Dyad(b)) 
LinearAlgebra.dot(a::PauliTypes{N}, b::FixedPhasePauli{N}) where N = get_coeff(a)*get_coeff(b) * dot(FixedPhasePauli(a), b) 
LinearAlgebra.dot(a::DyadTypes{N}, b::FixedPhasePauli{N}) where N = get_coeff(a)*get_coeff(b) * dot(Dyad(a), b) 