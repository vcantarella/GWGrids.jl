

"""
    get_linear_index(grid, l, r, c)

Returns the 1-based linear index for the cell at layer `l`, row `r`, column `c`.
Assumes standard Julia column-major ordering for a 3D array of size `(nlay, nrow, ncol)`.
"""
@inline function get_linear_index(grid::PlanarRegularGrid, l::Int, r::Int, c::Int)
    # Stride for row (moving 1 row adds nlay elements)
    # Stride for column (moving 1 column adds nlay*nrow elements)
    # Index = l + (r-1)*nlay + (c-1)*(nlay*nrow)
    return l + (r - 1) * grid.nlay + (c - 1) * (grid.nlay * grid.nrow)
end


"""
Converts a list of 3D (l, r, c) indices to a 1D linear index vector.
This function is GPU-aware.
"""
function _to_linear_indices(
    grid::PlanarRegularGrid, 
    locations::AbstractVector{<:Tuple{Int, Int, Int}}
)
    n_nodes = length(locations)
    
    # Create the output array (e.g., Vector{Int} or CuArray{Int})
    # We use grid.delr as a "template" to create an Int array on the correct device (CPU/GPU)
    indices_vec = similar(grid.delr, Int, n_nodes)
    
    # This loop is simple and will be slow on the CPU
    # but can be easily parallelized or run on GPU
    # For a real implementation, this would be a @kernel
    cpu_locations = locations 
    if !(locations isa Vector)
         cpu_locations = Array(locations) # Bring locations to CPU for iteration
    end
    
    indices_cpu = zeros(Int, n_nodes)
    for (i, (l, r, c)) in enumerate(cpu_locations)
        indices_cpu[i] = get_linear_index(grid, l, r, c)
    end
    
    # Copy the indices back to the original device
    # We use copyto! to ensure it goes to the right place if indices_vec is on GPU
    copyto!(indices_vec, indices_cpu)
    
    return indices_vec
end
"""
    compute layer, row and col from liner index, idx
"""
@inline function get_lrc(idx:Int,
    grid::PlanarRegularGrid)

    # Reverse the mapping: idx = (c - 1) * grid.nrow * grid.nlay + (r - 1) * grid.nlay + l
    idx0 = idx - 1
    
    l = (idx0 % grid.nlay) + 1
    
    rem_rc = idx0 ÷ grid.nlay
    r = (rem_rc % grid.nrow) + 1
    
    c = (rem_rc ÷ grid.nrow) + 1
    
    return (l, r, c)
end