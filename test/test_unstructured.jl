using Test
using GWGrids
using ExtendableGrids
using SimplexGridFactory

@testset "UnstructuredGrid Tests" begin
    # Create a 2D unstructured grid using ExtendableGrids
    # This creates a grid from a tensor product of 1D arrays
    X = [0.0, 1.0, 2.0]
    Y = [0.0, 1.0]
    egrid = simplexgrid(X, Y)
    
    # Wrap it in GWGrids.UnstructuredGrid
    ugrid = UnstructuredGrid(egrid)
    
    # Check grid_size. 3 * 2 = 6 nodes
    @test grid_size(ugrid) == (6,)
    
    # Check get_xyz for valid indices
    # Node 1 is (0.0, 0.0, 0.0)
    x, y, z = get_xyz(ugrid, 1)
    @test x ≈ 0.0
    @test y ≈ 0.0
    @test z ≈ 0.0
    
    # Node 2 is (1.0, 0.0, 0.0) based on ExtendableGrids tensor grid creation
    x, y, z = get_xyz(ugrid, CartesianIndex(2))
    @test x ≈ 1.0
    @test y ≈ 0.0
    @test z ≈ 0.0

    # Test FVM properties pre-computation
    @test ugrid.n_nodes == 6
    @test ugrid.n_cells > 0
    @test ugrid.n_faces > 0

    # Ensure arrays have correct lengths based on counts
    @test size(ugrid.coordinates, 2) == ugrid.n_nodes
    @test size(ugrid.cell_nodes, 2) == ugrid.n_cells
    @test size(ugrid.cell_centers, 2) == ugrid.n_cells
    @test length(ugrid.cell_volumes) == ugrid.n_cells
    
    @test size(ugrid.face_nodes, 2) == ugrid.n_faces
    @test size(ugrid.face_cells, 2) == ugrid.n_faces
    @test length(ugrid.face_areas) == ugrid.n_faces
    @test size(ugrid.face_normals, 2) == ugrid.n_faces
    @test size(ugrid.cell_faces, 2) == ugrid.n_cells
    
    # Verify that volumes are roughly correct (tensor grid is 2.0x1.0, total vol is 2.0)
    @test sum(ugrid.cell_volumes) ≈ 2.0
end
