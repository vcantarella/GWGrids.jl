using ExtendableGrids

grid = simplexgrid(0:1.0:2, 0:1.0:2)

for T in [FaceNodes, FaceCells, CellFaces, FaceGeometries, FaceVolumes, FaceNormals]
    try
        val = grid[T]
        println(T, " => ", typeof(val))
    catch e
        println(T, " => Not found/Error")
    end
end
