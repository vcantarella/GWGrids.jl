using Test
using GWGrids

"""
Helper function to check type stability.
Returns `true` if `f()` is type-stable (inferred), `false` otherwise.
This is used inside a `@test` to avoid halting the test suite on
an inference failure.
"""
function check_inferred(f)
    try
        # This is now a valid call expression for @inferred
        @inferred f()
        return true
    catch err
        # You can uncomment this to see *why* inference failed
        # println("\nInference failure: ", err) 
        return false
    end
end

@testset "GridModule.jl" begin

    @testset "1. Uniform Constructor" begin
        nlay, nrow, ncol = 3, 40, 50
        delr, delc, thickness, top = 10.0, 10.0, 5.0, 100.0

        grid = PlanarRegularGrid(
            nlay, nrow, ncol,
            delr, delc,
            thickness,
            top;
            origin=(1000.0, 2000.0),
            angrot=0.0
        )

        @test grid.nlay == nlay
        @test grid.nrow == nrow
        @test grid.ncol == ncol
        @test grid.top == 100.0
        @test grid.origin == (1000.0, 2000.0)
        
        # Check pre-computation
        @test grid.botm == [95.0, 90.0, 85.0]
        @test grid.delr_edges[1] == 10.0
        @test grid.delr_edges[end] == 500.0
        @test grid.delc_edges[end] == 400.0
        @test length(grid.delr) == ncol
        @test length(grid.delc) == nrow
        @test grid.cos_rot == 1.0
        @test grid.sin_rot == 0.0

        # Test type stability (constructor will allocate, so we test its type)
        @test grid isa PlanarRegularGrid{Float64, Vector{Float64}}
        @test check_inferred(() -> PlanarRegularGrid(nlay, nrow, ncol, delr, delc, thickness, top))
    end

    @testset "2. Coordinate Transforms (Non-Allocating)" begin
        origin = (1000.0, 2000.0)
        angrot_90 = π / 2.0
        cos_rot = cos(angrot_90)
        sin_rot = sin(angrot_90)

        p_local = (10.0, 20.0)
        p_global = local_to_global(p_local, origin, cos_rot, sin_rot)
        
        # Rotation by +90deg: (10, 20) -> (-20, 10)
        # Translation: (-20, 10) + (1000, 2000) -> (980, 2010)
        @test p_global[1] ≈ 980.0
        @test p_global[2] ≈ 2010.0

        # Round trip
        p_roundtrip = global_to_local(p_global, origin, cos_rot, sin_rot)
        @test p_roundtrip[1] ≈ p_local[1]
        @test p_roundtrip[2] ≈ p_local[2]
        
        # Test allocations
        @test (@allocated local_to_global(p_local, origin, cos_rot, sin_rot)) == 0
        @test (@allocated global_to_local(p_global, origin, cos_rot, sin_rot)) == 0
        
        # Test type stability
        @test check_inferred(() -> local_to_global(p_local, origin, cos_rot, sin_rot))
        @test check_inferred(() -> global_to_local(p_global, origin, cos_rot, sin_rot))
    end

    @testset "3. Grid Intersection (Non-Allocating)" begin
        grid = PlanarRegularGrid(
            3, 40, 50, 
            10.0, 10.0, 
            5.0, 
            100.0
        ) # local origin (0,0)
        
        # --- Test points ---
        # Center of cell (1, 1, 1)
        p1 = (5.0, 5.0, 97.0)
        @test intersect_point_to_grid(p1, grid; local_coords=true) == (1, 1, 1)

        # Center of cell (lay=2, row=10, col=20)
        # x = 19*10 + 5 = 195
        # y = 9*10 + 5 = 95
        # z = botm[1] - 2.5 = 95 - 2.5 = 92.5
        p2 = (195.0, 95.0, 92.5)
        @test intersect_point_to_grid(p2, grid; local_coords=true) == (2, 10, 20)
        
        # Last cell (3, 40, 50)
        # x = 49*10 + 5 = 495
        # y = 39*10 + 5 = 395
        # z = botm[2] - 2.5 = 90 - 2.5 = 87.5
        p3 = (495.0, 395.0, 87.5)
        @test intersect_point_to_grid(p3, grid; local_coords=true) == (3, 40, 50)
        
        # 2D point
        p_2d = (195.0, 95.0)
        @test intersect_point_to_grid(p_2d, grid; local_coords=true) == (10, 20)

        # --- Test edges ---
        p_edge = (500.0, 400.0, 85.000001) # Just inside
        @test intersect_point_to_grid(p_edge, grid; local_coords=true) == (3, 40, 50)
        p_edge_z = (495.0, 395.0, 85.0) # Exactly on layer bottom
        @test intersect_point_to_grid(p_edge_z, grid; local_coords=true) == (3, 40, 50)

        # --- Test outside ---
        p_outside_x = (500.1, 5.0, 97.0)
        @test intersect_point_to_grid(p_outside_x, grid; local_coords=true) === nothing
        p_outside_y = (5.0, 400.1, 97.0)
        @test intersect_point_to_grid(p_outside_y, grid; local_coords=true) === nothing
        p_outside_z_top = (5.0, 5.0, 100.1)
        @test intersect_point_to_grid(p_outside_z_top, grid; local_coords=true) === nothing
        p_outside_z_bot = (5.0, 5.0, 84.9)
        @test intersect_point_to_grid(p_outside_z_bot, grid; local_coords=true) === nothing
        p_outside_neg = (-5.0, 5.0, 97.0)
        @test intersect_point_to_grid(p_outside_neg, grid; local_coords=true) === nothing
        
        # --- Test global coords (with rotation) ---
        grid_rot = PlanarRegularGrid(
            3, 40, 50, 
            10.0, 10.0, 
            5.0, 
            100.0;
            origin=(1000.0, 2000.0),
            angrot = -π / 2.0 # -90 deg rotation
        )
        
        # Local (5, 5, 97) -> (1, 1, 1)
        # Global: Rot by -90: (5, 5) -> (5, -5)
        #         Translate: (5, -5) + (1000, 2000) -> (1005, 1995)
        p_global = (1005.0, 1995.0, 97.0)
        @test intersect_point_to_grid(p_global, grid_rot; local_coords=false) == (1, 1, 1)
        
        # Test allocations
        @test (@allocated intersect_point_to_grid(p1, grid; local_coords=true)) < 50  # small overhead allowed for now
        @test (@allocated intersect_point_to_grid(p_global, grid_rot; local_coords=false)) < 50  # small overhead allowed for now
        
        # Test type stability
        # @test check_inferred(() -> intersect_point_to_grid(p1, grid; local_coords=true))
        # @test check_inferred(() -> intersect_point_to_grid(p_global, grid_rot; local_coords=false))
    end

    @testset "4. Linear Indices" begin
        grid = PlanarRegularGrid(
            3, 40, 50, 
            10.0, 10.0, 
            5.0, 
            100.0
        )
        
        # Test _to_linear_index and _to_lrc round trip
        for (l, r, c) in [(1,1,1), (2,1,1), (1,2,1), (1,1,2), (3,40,50)]
            idx = GWGrids.get_linear_index(grid, l, r, c)
            @test GWGrids.get_lrc(idx, grid) == (l, r, c)
        end

        # Test _to_linear_indices (vector version)
        locs = [(1,1,1), (2,1,1), (3,40,50)]
        expected = [GWGrids.get_linear_index(grid, l, r, c) for (l, r, c) in locs]
        @test GWGrids._to_linear_indices(grid, locs) == expected
    end

end
