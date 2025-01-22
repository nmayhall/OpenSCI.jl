using OpenSCI
using Test

@testset "OpenSCI.jl" begin
    include("tests1.jl")
    include("test_paulimult.jl")
end
