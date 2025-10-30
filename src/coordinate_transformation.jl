function local_to_global(
    local_coords::Tuple{T, T},
    origin::Tuple{T, T},
    cos_rot::T,
    sin_rot::T
) where T
    # Unpack local coordinates
    x_local, y_local = local_coords
    x_origin, y_origin = origin

    # Apply rotation
    x_rotated = cos_rot * x_local - sin_rot * y_local
    y_rotated = sin_rot * x_local + cos_rot * y_local

    # Translate to global coordinates
    x_global = x_rotated + x_origin
    y_global = y_rotated + y_origin

    return (x_global, y_global)
end

"""
Transforms global (x, y) coordinates to local grid coordinates.
This function is non-allocating.
"""
function global_to_local(
    global_coords::Tuple{T, T},
    origin::Tuple{T, T},
    cos_rot::T,
    sin_rot::T
) where T
    # Unpack
    x_global, y_global = global_coords
    x_origin, y_origin = origin

    # Translate
    x_translated = x_global - x_origin
    y_translated = y_global - y_origin

    # Apply inverse rotation
    # cos(-a) = cos(a)
    # sin(-a) = -sin(a)
    x_local = cos_rot * x_translated + sin_rot * y_translated
    y_local = -sin_rot * x_translated + cos_rot * y_translated

    return (x_local, y_local)
end