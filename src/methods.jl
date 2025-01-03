using PauliOperators
using OpenSCI

function LinearAlgebra.rmul!(dρ::DyadSum{N,T}, L::Lindbladian{N}, ρ::DyadSum{N,T}, thresh=1e-16) where {N,T} 
    empty!(dρ)

    # Unitary part
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

    # Non-unitary
    for i in 1:length(γ)
        for (pauli,pcoeff) in L.L.ops
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
    end

    return
end

function Base.:*(L::Lindbladian{N}, ρ::DyadSum{N,T}) where {N,T}
    # out = DyadSum(N,T=T)
    # rmul!(out, L, ρ)
    dρ = L.H*ρ - ρ*L.H
    for i in 1:length(L.γ)
        Li = L.L[i]
        dρ += L.γ[i] * (Li * (ρ * Li'))
        LL = Li * Li'
        dρ -= 0.5*L.γ[i]*(LL*ρ + ρ*LL)
    end
    return dρ 
end
function Base.:*(L::Lindbladian{N}, ρd::Dyad{N}, T=ComplexF64) where N
    return L*DyadSum(ρd, T=T) 
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