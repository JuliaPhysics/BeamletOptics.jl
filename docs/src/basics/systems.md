# Optical systems

A collection of optical elements forms an optical system. Optical systems are used together with beams for the [`solve_system!`](@ref) function. The system types differ in how they store their objects, not in how they are solved: both run through the same tracing methods. A [`System`](@ref) keeps its objects and object groups in a vector, a [`StaticSystem`](@ref) in a tuple and flattens the groups when it is built. Currently the following types are available:

```@repl
using BeamletOptics # hide
BeamletOptics.list_subtypes(BeamletOptics.AbstractSystem); # hide
```

Refer to the **Tutorials** section for examples on how to define optical systems.

```@docs; canonical=false
System
StaticSystem
```

## Changing a system

The objects of a [`System`](@ref) can be changed after its construction, e.g. to build up a setup
step by step, starting from an empty system:

```julia
system = System()
push!(system, lens, mirror)   # objects or object groups
delete!(system, mirror)       # removes this object, compared by identity
pop!(system)                  # removes the last top-level object and returns it
popat!(system, 1)             # removes the first top-level object and returns it
```

A [`StaticSystem`](@ref) can not be changed.

```@docs; canonical=false
Base.push!(::System, ::Vararg{BeamletOptics.ObjectOrGroup})
Base.pop!(::System)
Base.popat!(::System, ::Integer)
Base.delete!(::System, ::BeamletOptics.ObjectOrGroup)
```

## Solving systems

Optical systems are solved non-sequentially by the [`solve_system!`](@ref) function: for every ray the solver searches the system for the object that is hit next. A component can shortcut this search for the following ray by returning a [`BeamletOptics.Hint`](@ref). This is explained in more detail in the section: [Tracing logic](@ref).

```@docs; canonical=false
solve_system!
```

Some objects store data during a solve, e.g. a [`Detector`](@ref) its hits. This data is kept between solves unless the system is initialized:

```@docs; canonical=false
initialize!
```

## Configuration

The numerical tolerances and defaults of the solver are global settings. Each has a getter without arguments, e.g. `get_sdf_miss_distance()`, whose docstring is listed below. The values are [Preferences.jl](https://github.com/JuliaPackaging/Preferences.jl) preferences that are read when BMO is loaded, so a new value takes effect after a restart of Julia.

| Getter | Default | Meaning |
| :--- | :---: | :--- |
| [`get_invariant_threshold`](@ref) | `1e-6` | limit in \[m\] of the complex invariant of an astigmatic beamlet |
| [`get_lagrange_tolerance`](@ref) | `0.5` | relative deviation of the Lagrange invariant of a beamlet |
| [`get_sdf_surface_threshold`](@ref) | `1e-9` | SDF value in \[m\] below which a ray counts as starting on the surface |
| [`get_sdf_raymarch_eps`](@ref) | `1e-10` | surface distance and minimum step of the ray marching in \[m\] |
| [`get_sdf_miss_distance`](@ref) | `1e15` | path length in \[m\] after which a marched ray is a miss |
| [`get_internal_reflection_threshold`](@ref) | `1e-6` | deviation of the reflectances from 1 below which a reflection is total, dimensionless |
| [`get_line_plane_intersection_threshold`](@ref) | `1e-6` | magnitude of the cosine between line and plane normal below which they are parallel, dimensionless |
| [`get_orthogonality_threshold`](@ref) | `1e-10` | magnitude of the dot product of unit vectors that still counts as orthogonal, e.g. for beam support axes, dimensionless |
| [`get_default_r_max`](@ref) | `100` | rays per beam, the `r_max` of [`solve_system!`](@ref) |
| [`get_default_depth_max`](@ref) | `100` | branching levels of the beam tree, the `depth_max` of [`solve_system!`](@ref) |
| [`get_power_cutoff`](@ref) | `0` | fraction of the power of the root beam below which a split-off child beam is not traced |
| [`get_default_wavelength`](@ref) | `1e-6` | wavelength in \[m\] |
| [`get_default_waist`](@ref) | `1e-3` | beam waist radius in \[m\] |
| [`get_default_power`](@ref) | `1e-3` | beam power in \[W\] |
| [`get_progress_threshold`](@ref) | `5.0` | runtime in \[s\] after which a progress bar is shown |

Only three of them have a setter: [`set_invariant_threshold!`](@ref), [`set_lagrange_tolerance!`](@ref) and [`set_progress_threshold!`](@ref). Any other preference is set by its key, which is the name of the getter without `get_`, via `BeamletOptics.Config.set_preference!`:

```julia
set_lagrange_tolerance!(0.2)
BeamletOptics.Config.set_preference!("sdf_miss_distance", 1e12)   # restart Julia afterwards
```

```@docs; canonical=false
get_invariant_threshold
set_invariant_threshold!
get_lagrange_tolerance
set_lagrange_tolerance!
get_sdf_surface_threshold
get_sdf_raymarch_eps
get_sdf_miss_distance
get_internal_reflection_threshold
get_line_plane_intersection_threshold
get_orthogonality_threshold
get_default_r_max
get_default_depth_max
get_power_cutoff
get_default_wavelength
get_default_waist
get_default_power
get_progress_threshold
set_progress_threshold!
BeamletOptics.Config.set_preference!
```