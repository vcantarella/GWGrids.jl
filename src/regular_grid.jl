struct PlanarRegularGrid{T<:AbstractFloat, V<:AbstractVector{T}}
    # --- Dimensions ---
    nlay::Int
    nrow::Int
    ncol::Int

    # --- Geometry ---
    delr::V
    delc::V
    top::T
    botm::V # Length NLAY, elevations of layer bottoms

    # --- Pre-computed for Lookups (The Fix) ---
    delr_edges::V  # = cumsum(delr)
    delc_edges::V  # = cumsum(delc)

    # --- Georeferencing ---
    origin::Tuple{T, T}
    angrot::T
    
    # --- Pre-computed for Rotation (The Fix) ---
    cos_rot::T
    sin_rot::T

    # Inner constructor
    function PlanarRegularGrid(
        delr::V, 
        delc::V, 
        top::T, 
        botm::V;
        origin::Tuple{T, T} = (zero(T), zero(T)),
        angrot::T = zero(T)
    ) where {T<:AbstractFloat, V<:AbstractVector{T}}

        nlay = length(botm)
        nrow = length(delc)
        ncol = length(delr)
        
        @assert nlay > 0 && nrow > 0 && ncol > 0 "Grid dimensions must be positive"

        # --- Pre-compute and store ---
        delr_edges = cumsum(delr)
        delc_edges = cumsum(delc)
        cos_rot = cos(angrot)
        sin_rot = sin(angrot)
        # ---
        
        # Note: If `V` is a CuArray, cumsum() will return a CuArray. This is correct!
        new{T, V}(
            nlay, nrow, ncol, 
            delr, delc, top, botm,
            delr_edges, delc_edges,
            origin, angrot,
            cos_rot, sin_rot
        )
    end
end

"""
Constructs a PlanarRegularGrid with uniform cell spacing.

This is a convenience constructor for creating a grid where all cells have
the same length, width, and thickness.

# Arguments
- `nlay::Int`: Number of layers.
- `nrow::Int`: Number of rows.
- `ncol::Int`: Number of columns.
- `Δx::Real`: Uniform column spacing (size in x-direction, Δx).
- `Δy::Real`: Uniform row spacing (size in y-direction, Δy).
- `Δz::Real`: Uniform thickness for all layers (size in z-direction, Δz).
- `top::Real`: The elevation of the grid top (a single planar value).

# Keyword Arguments
- `origin::Tuple{<:Real, <:Real} = (0.0, 0.0)`: Global (x, y) coordinates of the grid's (0, 0) local corner.
- `angrot::Real = 0.0`: Counter-clockwise rotation angle in radians.
"""
function PlanarRegularGrid(
    nlay::Int, nrow::Int, ncol::Int;
    delr::Real,                 # Uniform column spacing (delta-x)
    delc::Real,                 # Uniform row spacing (delta-y)
    layer_thickness::Real,      # Uniform layer thickness (delta-z)
    top::Real,                  # Top elevation of the grid
    origin::Tuple{<:Real, <:Real} = (0.0, 0.0),
    angrot::Real = 0.0
)
    # Determine the promotion type for all floating-point geometry
    T = promote_type(
        typeof(delr), typeof(delc), typeof(layer_thickness), typeof(top), 
        eltype(origin), typeof(angrot)
    )
    
    # Ensure it's a float type, defaulting to Float64
    T_float = T <: AbstractFloat ? T : Float64

    # Convert all inputs to the promoted float type
    delr_t = T_float(delr)
    delc_t = T_float(delc)
    layer_thickness_t = T_float(layer_thickness)
    top_t = T_float(top)
    origin_tup = (T_float(origin[1]), T_float(origin[2]))
    angrot_val = T_float(angrot)

    # --- Create the geometry arrays ---
    
    # Fill standard CPU `Vector`s with the uniform values
    delr_vec = fill(delr_t, ncol)
    delc_vec = fill(delc_t, nrow)
    
    # Calculate planar bottom elevations
    # botm[1] = top - 1*thickness
    # botm[2] = top - 2*thickness
    # ...
    botm_vec = [top_t - k * layer_thickness_t for k in 1:nlay]

    # --- Call the main (inner) constructor ---
    # This will create a grid using standard `Vector`s for CPU usage.
    # The inner constructor will automatically pre-compute the
    # `delr_edges`, `delc_edges`, `cos_rot`, and `sin_rot` fields.
    return PlanarRegularGrid(
        delr_vec, 
        delc_vec, 
        top_t, 
        botm_vec; 
        origin=origin_tup, 
        angrot=angrot_val
    )
end


"""
Holds the flow solution (head and face flows) on a grid.

The array types `ArrT3D` can be `Array` (for CPU) or `CuArray` (for GPU).
"""
struct FlowSolution{
    T<:AbstractFloat, 
    ArrT3D<:AbstractArray{T, 3}, 
    G<:PlanarRegularGrid{T}
}
    grid::G
    head::ArrT3D  # Cell-centered head, size (nlay, nrow, ncol)

    # Face flows stored in a NamedTuple for clarity
    flows::NamedTuple{(:right, :front, :lower), NTuple{3, ArrT3D}}

    # Constructor
    function FlowSolution(
        grid::G, 
        head::ArrT3D, 
        # This constructor signature is flexible and correct
        flows::NamedTuple{(:right, :front, :lower), <:NTuple{3, ArrT3D}}
    ) where {T<:AbstractFloat, ArrT3D<:AbstractArray{T, 3}, G<:PlanarRegularGrid{T}}
        
        dims = (grid.nlay, grid.nrow, grid.ncol)
        @assert size(head) == dims "Head array has incorrect dimensions"
        @assert size(flows.right) == dims "Flow-Right-Face has incorrect dimensions"
        @assert size(flows.front) == dims "Flow-Front-Face has incorrect dimensions"
        @assert size(flows.lower) == dims "Flow-Lower-Face has incorrect dimensions"
        
        new{T, ArrT3D, G}(grid, head, flows)
    end
end