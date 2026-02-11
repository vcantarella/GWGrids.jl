module GWGrids
    include("regular_grid.jl")
    include("coordinate_transformation.jl")
    include("intersects.jl")
    include("linear_indexes.jl")
    include("grid_types.jl")
    export PlanarRegularGrid, intersect_point_to_grid, global_to_local, local_to_global
    export get_linear_index, get_lrc
    export AbstractGWGrid
    export UniformGrid, RectilinearGrid, StructuredGrid, UnstructuredGrid
    export get_xyz, grid_size

end # module GWGrids
