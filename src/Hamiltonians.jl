using PauliOperators
using OpenSCI

function heisenberg_1D(N, Jx, Jy, Jz)
    H = PauliSum(N, ComplexF64)
    for i in 0:N-1
        H += -2*Jx * Pauli(N, X=[i+1,(i+1)%(N)+1])
        H += -2*Jy * Pauli(N, Y=[i+1,(i+1)%(N)+1])
        H += -2*Jz * Pauli(N, Z=[i+1,(i+1)%(N)+1])
    end 
    return H
end 