struct PlanarRegularGrid{T<:AbstractFloat, V<:AbstractVector{T}, I<:AbstractVector{Int}}
    # --- Dimensions ---
    nlay::Int
    nrow::Int
    ncol::Int

    # --- Geometry ---
    delr::V
    delc::V
    top::T
    botm::V # Length NLAY, elevations of layer bottoms
    
    # --- Properties ---
    li::I # Layer Type Indicator (e.g., 0=confined, 1=convertible)

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
        li::Union{Nothing, I} = nothing,
        origin::Tuple{T, T} = (zero(T), zero(T)),
        angrot::T = zero(T)
    ) where {T<:AbstractFloat, V<:AbstractVector{T}, I<:AbstractVector{Int}}

        nlay = length(botm)
        nrow = length(delc)
        ncol = length(delr)
        
        @assert nlay > 0 && nrow > 0 && ncol > 0 "Grid dimensions must be positive"
        
        # Default li to zeros (confined) if not provided
        # We need to ensure the vector type matches the device of V if possible,
        # but usually li is just a small 1D vector. 
        # For simplicity, if V is a GPU array, the user should pass a GPU Int array for li.
        # If defaulted, we'll make a CPU Vector{Int}.
        
        local_li::I = if li === nothing
             # If V is a CuArray, we can't easily infer the corresponding Int CuArray type 
             # without CUDA.jl dependencies here.
             # Ideally, the user provides it.
             # Fallback: assume standard Vector{Int} and let the user convert later if needed,
             # OR strictly require it if V is not a standard Vector.
             if V <: Vector
                 zeros(Int, nlay)
             else
                 # If we are on GPU, creating a CPU array might be mixed usage.
                 # But let's create a CPU one and let the user handle transfer if needed,
                 # or error out if strictness is preferred.
                 # For now: safe default.
                 zeros(Int, nlay)
             end
        else
            li
        end
        
        @assert length(local_li) == nlay "Length of li must match nlay"

        # --- Pre-compute and store ---
        delr_edges = cumsum(delr)
        delc_edges = cumsum(delc)
        cos_rot = cos(angrot)
        sin_rot = sin(angrot)
        # ---
        
        new{T, V, typeof(local_li)}(
            nlay, nrow, ncol, 
            delr, delc, top, botm,
            local_li,
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
- `li::Vector{Int}`: Layer indicator array (length nlay). Defaults to zeros.
- `origin::Tuple{<:Real, <:Real} = (0.0, 0.0)`: Global (x, y) coordinates of the grid's (0, 0) local corner.
- `angrot::Real = 0.0`: Counter-clockwise rotation angle in radians.
"""
function PlanarRegularGrid(
    nlay::Int, nrow::Int, ncol::Int,
    Δx::Real,                 # Uniform column spacing (delta-x)
    Δy::Real,                 # Uniform row spacing (delta-y)
    Δz::Real,                 # Uniform thickness for all layers (delta-z)
    top::Real;                # Top elevation of the grid
    li::Union{Vector{Int}, Nothing} = nothing,
    origin::Tuple{<:Real, <:Real} = (0.0, 0.0),
    angrot::Real = 0.0
)
    # Determine the promotion type for all floating-point geometry
    T = promote_type(
        typeof(Δx), typeof(Δy), typeof(Δz), typeof(top), 
        eltype(origin), typeof(angrot)
    )
    
    # Ensure it's a float type, defaulting to Float64
    T_float = T <: AbstractFloat ? T : Float64

    # Convert all inputs to the promoted float type
    Δx_t = T_float(Δx)
    Δy_t = T_float(Δy)
    Δz_t = T_float(Δz)
    top_t = T_float(top)
    origin_tup = (T_float(origin[1]), T_float(origin[2]))
    angrot_val = T_float(angrot)

    # --- Create the geometry arrays ---
    
    # Fill standard CPU `Vector`s with the uniform values
    delr_vec = fill(Δx_t, ncol)
    delc_vec = fill(Δy_t, nrow)

    # Calculate planar bottom elevations
    botm_vec = [top_t - k * Δz_t for k in 1:nlay]
    
    # Default li if not provided
    li_vec = (li === nothing) ? zeros(Int, nlay) : li

    # --- Call the main (inner) constructor ---
    return PlanarRegularGrid(
        delr_vec, 
        delc_vec, 
        top_t, 
        botm_vec; 
        li=li_vec,
        origin=origin_tup, 
        angrot=angrot_val
    )
end



#     T<:AbstractFloat, 
#     ArrT3D<:AbstractArray{T, 3}, 
#     G<:PlanarRegularGrid{T}
# }
#     grid::G
#     head::ArrT3D  # Cell-centered head, size (nlay, nrow, ncol)

#     # Face flows stored in a NamedTuple for clarity
#     flows::NamedTuple{(:right, :front, :lower), NTuple{3, ArrT3D}}

#     # Constructor
#     function FlowSolution(
#         grid::G, 
#         head::ArrT3D, 
#         # This constructor signature is flexible and correct
#         flows::NamedTuple{(:right, :front, :lower), <:NTuple{3, ArrT3D}}
#     ) where {T<:AbstractFloat, ArrT3D<:AbstractArray{T, 3}, G<:PlanarRegularGrid{T}}
        
#         dims = (grid.nlay, grid.nrow, grid.ncol)
#         @assert size(head) == dims "Head array has incorrect dimensions"
#         @assert size(flows.right) == dims "Flow-Right-Face has incorrect dimensions"
#         @assert size(flows.front) == dims "Flow-Front-Face has incorrect dimensions"
#         @assert size(flows.lower) == dims "Flow-Lower-Face has incorrect dimensions"
        
#         new{T, ArrT3D, G}(grid, head, flows)
#     end
# end