using PauliOperators
using OpenSCI

function LinearAlgebra.rmul!(dρ::DyadSum{N,T}, L::Lindbladian{N}, ρ::DyadSum{N,T}, thresh=1e-16) where {N,T} 
    empty!(dρ)
    for (pauli,pcoeff) in L.H.ops
       
        for (dyad,dcoeff) in ρ
            sd = pauli * dyad
            coeff = sd.coeff * pcoeff * dcoeff
            abs2(coeff) > thresh || continue
            sum!(dρ, coeff*sd.dyad)
            
            sd = dyad * pauli
            coeff = sd.coeff * pcoeff * dcoeff
            abs2(coeff) > thresh || continue
            sum!(dρ, -coeff*sd.dyad)
        end
    end
    return
end

function Base.:*(L::Lindbladian{N}, ρ::DyadSum{N,T}) where {N,T}
    out = DyadSum(N,T=T)
    rmul!(out, L, ρ)
    return out
end
function Base.:*(L::Lindbladian{N}, ρ::Dyad{N}, T=ComplexF64) where N
    out = DyadSum(N,T=T)
    rmul!(out, L, DyadSum(ρ, T=T))
    return out
end

Base.vec(d::Dyad{N}) where N = 1 + d.ket.v + d.bra.v*(BigInt(2)^N)

function Base.Matrix(L::Lindbladian{N}) where N
    Lmat = zeros(ComplexF64, 4^N, 4^N)
    for kr in 0:2^N-1
        for br in 0:2^N-1
            dyad_r = Dyad(N, kr, br)
            
            # rmul!(σ, L, DyadSum(dyad_r))
            σ = L * dyad_r        
            for (dyad_l, coeff) in σ
                Lmat[vec(dyad_l), vec(dyad_r)] = coeff
            end
        end
    end
    return Lmat
end