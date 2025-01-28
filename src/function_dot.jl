# function LinearAlgebra.dot(a::PauliBasis{N}, b::DyadBasis{N}) where N
#     i1 = iseven(count_ones(a.z & b.bra.v)) 
#     i2 = b.bra.v == (a.x ⊻ b.ket.v)
#     return -(-1)^i1*i2 
# end

# LinearAlgebra.dot(a::PauliBasis{N}, b::PauliBasis{N}) where N = a == b
# LinearAlgebra.dot(a::Dyad{N}, b::Dyad{N}) where N = (a.ket == b.bra) && (a.bra == b.ket)

# PauliTypes{N} = Union{Pauli{N}, ScaledPauli{N}}
# DyadTypes{N} = Union{Dyad{N}, ScaledDyad{N}}

# PauliOperators.get_coeff(p::PauliBasis) = 1
# PauliOperators.get_coeff(p::Dyad) = 1
# PauliOperators.get_coeff(p::ScaledDyad) = p.coeff
# # PauliOperators.PauliBasis(p::PauliBasis) = p 
# PauliOperators.PauliBasis(p::Pauli) = p.pauli 
# PauliOperators.PauliBasis(p::ScaledPauli) = p.pauli 
# # PauliOperators.Dyad(d::Dyad) = d
# PauliOperators.Dyad(d::ScaledDyad) = d.dyad

# LinearAlgebra.dot(a::PauliTypes{N}, b::DyadTypes{N}) where N = get_coeff(a)*get_coeff(b) * dot(PauliBasis(a), Dyad(b)) 
# LinearAlgebra.dot(a::DyadTypes{N}, b::DyadTypes{N}) where N = get_coeff(a)*get_coeff(b) * dot(Dyad(a), Dyad(b)) 
# LinearAlgebra.dot(a::DyadTypes{N}, b::PauliTypes{N}) where N = get_coeff(a)*get_coeff(b) * dot(Dyad(a), PauliBasis(b)) 
# LinearAlgebra.dot(a::PauliTypes{N}, b::PauliTypes{N}) where N = get_coeff(a)*get_coeff(b) * dot(PauliBasis(a), PauliBasis(b)) 
# LinearAlgebra.dot(a::PauliBasis{N}, b::PauliTypes{N}) where N = get_coeff(a)*get_coeff(b) * dot(a, PauliBasis(b)) 
# LinearAlgebra.dot(a::PauliBasis{N}, b::DyadTypes{N}) where N = get_coeff(a)*get_coeff(b) * dot(a, Dyad(b)) 
# LinearAlgebra.dot(a::PauliTypes{N}, b::PauliBasis{N}) where N = get_coeff(a)*get_coeff(b) * dot(PauliBasis(a), b) 
# LinearAlgebra.dot(a::DyadTypes{N}, b::PauliBasis{N}) where N = get_coeff(a)*get_coeff(b) * dot(Dyad(a), b) 