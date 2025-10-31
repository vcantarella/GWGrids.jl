
"""
Converts a 3D (layer, row, col) index to a 1D linear index.
"""
function _to_linear_index(grid::PlanarRegularGrid, l::Int, r::Int, c::Int)
    # Perform bounds checking
    if !(1 <= l <= grid.nlay && 1 <= r <= grid.nrow && 1 <= c <= grid.ncol)
        error("Grid index ($l, $r, $c) is out of bounds for grid size ($grid.nlay, $grid.nrow, $grid.ncol)")
    end
    
    # Calculate linear index (column-major order is default in Julia)
    # But MODFLOW uses (lay, row, col) which is more like (dim3, dim2, dim1)
    # We'll assume a layout of (nlay, nrow, ncol)
    # idx = c + (r-1)*grid.ncol + (l-1)*grid.ncol*grid.nrow
    
    # Let's use Julia's built-in (and optimized) function for this.
    # Note: `CartesianIndex` is (i, j, k) which maps to (row, col, lay)
    # It's often safer to define the mapping explicitly.
    # Assuming (nlay, nrow, ncol) storage:
    idx = (c-1) * grid.nrow * grid.nlay + (r-1) * grid.nlay + l
    
    # Let's double-check. (1,1,1) -> 1. (2,1,1) -> 2. (1,2,1) -> nlay + 1.
    # This seems like a (lay, row, col) mapping.
    
    # A safer, more standard Julia approach (assumes (row, col, lay) storage)
    # idx = (l-1) * grid.nrow * grid.ncol + (c-1) * grid.nrow + r
    
    # Let's stick to the MODFLOW (lay, row, col) layout, but stored
    # in a Julia (nlay, nrow, ncol) array.
    # LinearIndex = (c-1)*nrow*nlay + (r-1)*nlay + l
    # This is for (nlay, nrow, ncol)
    # Let's test:
    # (1,1,1) -> 0*... + 0*... + 1 = 1
    # (2,1,1) -> 0*... + 0*... + 2 = 2
    # (1,2,1) -> 0*... + 1*nlay + 1 = nlay + 1
    # (1,1,2) -> 1*nrow*nlay + 0*nlay + 1 = nrow*nlay + 1
    # This seems correct for an array A[l, r, c]
    
    linear_idx = (c - 1) * grid.nrow * grid.nlay + (r - 1) * grid.nlay + l
    return linear_idx
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
    
    # Get the correct array type (Vector or CuArray) from the grid
    # We use grid.delr as a "template" for the array type
    ArrayType = typeof(grid.delr)
    
    # Create the output array (e.g., Vector{Int} or CuArray{Int})
    # We must use `similar` to create an Int array on the correct device (CPU/GPU)
    indices_vec = similar(ArrayType, Int, n_nodes)
    
    # This loop is simple and will be slow on the CPU
    # but can be easily parallelized or run on GPU
    # For a real implementation, this would be a @kernel
    cpu_locations = locations #
    if !(locations isa Vector)
         cpu_locations = Array(locations) # Bring locations to CPU for iteration
    end
    
    indices_cpu = zeros(Int, n_nodes)
    for (i, (l, r, c)) in enumerate(cpu_locations)
        indices_cpu[i] = _to_linear_index(grid, l, r, c)
    end
    
    # Copy the indices back to the original device
    indices_vec = ArrayType(indices_cpu)
    
    return indices_vec
end