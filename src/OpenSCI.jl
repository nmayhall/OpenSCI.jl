module OpenSCI

using PauliOperators
using LinearAlgebra
using Printf

include("types.jl")
include("methods.jl")
include("sci.jl")
include("Hamiltonians.jl")

export Lindbladian
export SubspaceDissipator 
export SparseDyadVectors
export todense
export selected_ci
export multiply
export build_subspace_L
export add_hamiltonian! 
export add_channel_amplitude_damping! 
export add_channel_dephasing! 
export add_channel_depolarizing! 

end
