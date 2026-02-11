@testset "New Grid Types" begin

    @testset "UniformGrid" begin
        # 1. Initialization
        origin = (0.0, 0.0, 0.0)
        spacing = (1.0, 2.0, 0.5)
        dims = (10, 5, 20)
        
        grid = UniformGrid(origin, spacing, dims)
        
        # 2. Test grid_size
        @test grid_size(grid) == (10, 5, 20)
        
        # 3. Test get_xyz (Math based)
        # First cell (origin)
        I_start = CartesianIndex(1, 1, 1)
        x, y, z = get_xyz(grid, I_start)
        @test (x, y, z) == (0.0, 0.0, 0.0)
        
        # Arbitrary cell: (2, 2, 2) -> (0 + 1*1, 0 + 1*2, 0 + 1*0.5) = (1.0, 2.0, 0.5)
        I_mid = CartesianIndex(2, 2, 2)
        x, y, z = get_xyz(grid, I_mid)
        @test (x, y, z) == (1.0, 2.0, 0.5)
    end

    @testset "RectilinearGrid" begin
        # 1. Initialization
        x_coords = [0.0, 1.0, 5.0, 10.0]
        y_coords = [0.0, 2.0, 4.0]
        z_coords = [10.0, 20.0]
        
        grid = RectilinearGrid(x_coords, y_coords, z_coords)
        
        # 2. Test grid_size
        @test grid_size(grid) == (4, 3, 2)
        
        # 3. Test get_xyz (Direct lookup)
        # Look up index (3, 2, 1) -> x[3], y[2], z[1]
        I = CartesianIndex(3, 2, 1)
        xi, yi, zi = get_xyz(grid, I)
        
        @test xi == 5.0
        @test yi == 2.0
        @test zi == 10.0
    end

    @testset "StructuredGrid" begin
        # 1. Initialization
        # Regular X and Y
        x_coords = [10.0, 20.0]
        y_coords = [5.0, 15.0]
        
        # Deformed Z (2x2x2 array)
        z_array = zeros(2, 2, 2)
        z_array[1, 1, 1] = 50.0
        z_array[2, 2, 1] = 99.0 # Specific point to test
        
        grid = StructuredGrid(x_coords, y_coords, z_array)
        
        # 2. Test grid_size (Should match Z array size)
        @test grid_size(grid) == (2, 2, 2)
        
        # 3. Test get_xyz (Hybrid lookup)
        # Case A: (1, 1, 1) -> x[1], y[1], z[1,1,1]
        xa, ya, za = get_xyz(grid, CartesianIndex(1, 1, 1))
        @test (xa, ya, za) == (10.0, 5.0, 50.0)
        
        # Case B: (2, 2, 1) -> x[2], y[2], z[2,2,1]
        xb, yb, zb = get_xyz(grid, CartesianIndex(2, 2, 1))
        @test (xb, yb, zb) == (20.0, 15.0, 99.0)
    end

    @testset "UnstructuredGrid" begin
        # 1. Initialization (List of points)
        x_pts = [1.0, 2.0, 3.0]
        y_pts = [10.0, 20.0, 30.0]
        z_pts = [100.0, 200.0, 300.0]
        
        grid = UnstructuredGrid(x_pts, y_pts, z_pts)
        
        # 2. Test grid_size
        @test grid_size(grid) == (3,)
        
        # 3. Test get_xyz (Linear lookup)
        # Case A: Integer Index (common in 1D kernels)
        xi, yi, zi = get_xyz(grid, 2)
        @test (xi, yi, zi) == (2.0, 20.0, 200.0)
        
        # Case B: CartesianIndex{1} (common in KA.jl kernels)
        xi_c, yi_c, zi_c = get_xyz(grid, CartesianIndex(3))
        @test (xi_c, yi_c, zi_c) == (3.0, 30.0, 300.0)
    end

end