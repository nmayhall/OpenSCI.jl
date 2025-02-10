using OrderedCollections


SparseDyadVectors{N,T} = OrderedDict{DyadBasis{N}, Vector{T}} 

function SparseDyadVectors(ds::DyadSum{N,T}) where {N,T}
    sdv = OrderedDict{DyadBasis{N}, Vector{T}}()
    for (dyad,coeff) in ds
        sdv[dyad] = [coeff] 
    end
    return sdv
end

function clip!(ps::SparseDyadVectors{N,T}; thresh=1e-16) where {N,T}
    filter!(p->maximum(abs.(p.second)) > thresh, ps.ops)
end

Base.Vector(sdv::SparseDyadVectors{N,T}; state=1)  where {N,T} = [i[state] for i in values(sdv)]
function Base.Matrix(sdv::SparseDyadVectors{N,T})  where {N,T}
    ni, nj = size(sdv)
    out = zeros(T, ni, nj)
    i = 1
    for (d,coeffs) in sdv
        for j in 1:nj 
            out[i,j] = coeffs[j]
        end
        i += 1
    end
    return out
end

function Base.size(sdv::SparseDyadVectors)
    ni = length(sdv)
    for (d,c) in sdv
        return (ni, length(c)) 
    end
end

function Base.sum!(sdv::SparseDyadVectors{N,T}, dyad::DyadBasis{N}, coeffs::Vector{T}) where {N,T}
    if haskey(sdv, dyad)
        sdv[dyad] .+= coeffs
    else
        sdv[dyad] = coeffs
    end
end

function todense(sdv::SparseDyadVectors{N,T}) where {N,T}
    out = zeros(T, 4^N, size(sdv)[2])
    for (d,c) in sdv
        out[index(d),:] .= c
    end
    return out
end