using Test
using GWGrids
using JET

@testset "GWGrids.jl" begin
    include("test_grids.jl")
    include("test_new_grids.jl")
    include("test_unstructured.jl")
end

@testset "JET tests" begin
    JET.test_package(GWGrids; target_modules=(GWGrids,))
end