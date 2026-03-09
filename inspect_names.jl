using ExtendableGrids
println(filter(x -> occursin("Cell", String(x)) || occursin("Face", String(x)), names(ExtendableGrids)))
