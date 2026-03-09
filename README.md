# GWGrids.jl
minimal grids for finite differences groundwater flow models (raw)

## Features

- **Regular Grids**: `UniformGrid`, `RectilinearGrid`, and `StructuredGrid`.
- **Unstructured Grids**: `UnstructuredGrid` is natively supported via integration with `ExtendableGrids.jl` and `SimplexGridFactory.jl`. Allowing flexible unstructured mesh creations directly from these established finite element and finite volume packages. 

## Testing and Quality Assurance

- Extensive tests are implemented via Julia's `Test` module.
- Type stabilities and potential issues are thoroughly checked with `JET.jl`.
