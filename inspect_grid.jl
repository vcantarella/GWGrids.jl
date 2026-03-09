using ExtendableGrids

grid = simplexgrid(0:1.0:2, 0:1.0:2)
println("Keys:")
for k in keys(grid)
    println("  ", k, " => ", typeof(grid[k]))
end

# Check for specific items
println("\nSpecific components:")
for T in [CellGeometries, CellNodes, CellVolumes, BFaceNodes, BFaceGeometries, BFaceRegions, 
          Coordinates, CellRegions, BFaceVolumes]
    try
        val = grid[T]
        println("  ", T, " => ", typeof(val))
    catch e
        println("  ", T, " => Not found/Error")
    end
end
