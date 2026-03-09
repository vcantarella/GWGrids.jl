using ExtendableGrids

grid = simplexgrid(0:1.0:2, 0:1.0:2)

for T in [CellCenters, CellMidpoints, CellBarycenters]
    try
        val = grid[T]
        println(T, " => ", typeof(val))
    catch e
        println(T, " => Not found/Error")
    end
end
