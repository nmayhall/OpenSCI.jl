
"""
    index(k::Ket)

Return location of this basis vector in the basis of `KetBitStrings`
"""
index(k::Ket) = 1 + k.v

"""
    index(d::DyadBasis{N}) where N

Return location of this basis vector in the basis of `Dyad`'s
"""
index(d::DyadBasis{N}) where N = 1 + d.ket.v + d.bra.v*(2^N)
index(d::Dyad) = index(DyadBasis(d)) 

"""
    index(p::PauliBasis{N}) where N

Return location of this basis vector in the basis of `FixedPhasePauli`'s
"""
index(p::PauliBasis{N}) where N = 1 + p.z + p.x*(2^N)
index(d::Pauli) = index(PauliBasis(d)) 