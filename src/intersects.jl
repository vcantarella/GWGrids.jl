"""
Finds the (layer, row, col) index for a 3D point.
returns nothing if the point is outside the grid.
If the point is 2D (no z-coordinate), returns (row, col).
"""
function intersect_point_to_grid(
    point::P, # Can be a Tuple, Vector, SVector, etc.
    grid::PlanarRegularGrid{T};
    local_coords::Bool = false
) where {T, P}
    
    # --- 1. Get Local Coordinates (Non-mutating) ---
    x_global, y_global = point[1], point[2]

    x, y = if local_coords
        (x_global, y_global)
    else
        global_to_local(
            (x_global, y_global), 
            grid.origin, 
            grid.cos_rot, 
            grid.sin_rot
        )
    end
    
    # --- 2. Find Column (X) and Row (Y) ---
    
    # Check boundaries (x < 0 or y < 0)
    if x < zero(T) || y < zero(T)
        return nothing
    end

    # Find column (binary search)
    col = searchsortedfirst(grid.delr_edges, x)
    if col > grid.ncol
        return nothing # Point is outside grid in x-direction
    end

    # Find row (binary search)
    row = searchsortedfirst(grid.delc_edges, y)
    if row > grid.nrow
        return nothing # Point is outside grid in y-direction
    end

    # --- 3. Find Layer (Z) ---
    if length(point) < 3
        return (row, col) # Return 2D index
    end

    z = point[3]
    if z > grid.top
        return nothing # Point is above the grid
    end
    
    # Find layer (binary search, in reverse order)
    # This finds the first index `k` where botm[k] < z
    layer = searchsortedfirst(grid.botm, z, rev=true)
    
    if layer > grid.nlay
        return nothing # Point is below the last layer
    end

    return (layer, row, col)
end