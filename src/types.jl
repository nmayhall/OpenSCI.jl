using PauliOperators
using OrderedCollections
# using LinearAlgebra
# using Printf
# using Dictionaries
# using StaticArrays





# function Base.display(sdv::DyadSum{N,T})  where {N,T}
#     idx = 1
#     for (d,coeffs) in sdv
#         @printf(" %6i %12s", idx, d)
#         for coeff in coeffs
#             @printf(" % 10.8f % 10.8fi", real(coeff), imag(coeff))
#         end
#         @printf("\n")
#         idx += 1
#     end
# end
