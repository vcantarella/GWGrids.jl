module GWGrids
    include("regular_grid.jl")
    include("coordinate_transformation.jl")
    include("intersects.jl")
    include("linear_indexes.jl")
    export PlanarRegularGrid, intersect_point_to_grid, global_to_local, local_to_global
    export get_linear_index, get_lrc
end # module GWGrids
