module GWGrids
    include("regular_grid.jl")
    include("coordinate_transformation.jl")
    include("intersects.jl")
    export PlanarRegularGrid, intersect_point_to_grid, global_to_local, local_to_global
end # module GWGrids
