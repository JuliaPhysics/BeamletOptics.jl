```@setup rays
dir = joinpath(@__DIR__, "..", "assets")

Main.DocUtils.conditional_include(joinpath(dir, "meshes.jl"))
```

# Meshes

For flat surfaces and placeholder objects, BMO offers a mesh-based geometry representation via the `Mesh` type. Meshes offer an easy but inaccurate way to represent geometry via a tessellation based surface interpolation with triangles. In general it is not recommended to use meshes for anything other than flat surfaces or `.stl` imports for visualization purposes. The point of intersection is calculated via the [Moeller Trumbore algorithm](@ref) [Moeller:2005](@cite).

Meshes can be loaded from `.stl` files via the `FileIO` package.

```@docs; canonical=false
BeamletOptics.Mesh
BeamletOptics.Mesh(::Any)
```

## Moeller Trumbore algorithm

The Moeller-Trumbore-algorithm tests the intersection of a ray with the mesh by testing against each individual face.

```@docs; canonical=false
BeamletOptics.MoellerTrumboreAlgorithm
```

Below a visualization is provided for a Benchy mesh. The ray, colored in blue, intersects multiple faces of the mesh, colored in red. The shortest face intersection is assumed to be the correct one. Note that the normal vector at the point of intersection is calculated for the entire face. This makes the mesh-based geometry representation unsuitable for curved geometries, e.g. lens surfaces, since refraction and reflection calculations require accurate normal vectors.

![Moeller Trumbore algorithm showcase](mtalgorithm.png)