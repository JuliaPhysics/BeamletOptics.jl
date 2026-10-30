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

The test for rays parallel to the face is relative to the size of the face. The determinant `Det = E1 ⋅ (d × E2)` of the edges `E1`, `E2` and the ray direction `d` is a triple product: it scales with the area of the face, `|Det| = |E1| |E2| sin α |cos θ|` with the corner angle `α` between the edges and the angle `θ` between `d` and the face normal. The ray is skipped if `|Det| ≤ kϵ |E1| |E2|`, which does not depend on the edge lengths. A degenerate face, with `Det = 0`, is skipped as well.

The vertices of a face, here and in `normal3d(mesh, fID)`, are converted to points in the number type of the ray or the mesh and type-asserted, since their type cannot be inferred for every mesh type. For the same reason the distance `t` of a hit is declared in the number type of the ray. This keeps both functions type-stable. [`BeamletOptics.CircularFlatMesh`](@ref) stores its vertices in the number type of the radius, e.g. dual numbers in automatic differentiation, so that the mesh can be differentiated with respect to the radius.

Below a visualization is provided for a Benchy mesh. The ray, colored in blue, intersects multiple faces of the mesh, colored in red. The shortest face intersection is assumed to be the correct one. Note that the normal vector at the point of intersection is calculated for the entire face. This makes the mesh-based geometry representation unsuitable for curved geometries, e.g. lens surfaces, since refraction and reflection calculations require accurate normal vectors.

![Moeller Trumbore algorithm showcase](mtalgorithm.png)