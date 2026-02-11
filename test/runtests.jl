using Test
using GWGrids

@testset "GWGrids.jl" begin
    include("test_grids.jl")
    include("test_new_grids.jl")
end