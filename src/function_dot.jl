
MyOps{N,T} = Union{Dyad{N}, DyadSum{N,T}, Pauli{N}, PauliSum{N,T}}

LinearAlgebra.dot(a::MyOps, b::MyOps) = tr(a*b)
