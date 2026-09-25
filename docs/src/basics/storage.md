# Saving and loading setups

A complete setup, i.e. one or more [`System`](@ref)s together with the beam sources that
illuminate them, paired the same way as for the GUI's `live_view(system => source, ...)` or as
`system => (source1, source2)` for a system traced by several sources, can be written to and read back from a single archive file (`*.bmo` by convention). Objects
and sources that occur several times, e.g. a source that illuminates two systems, are
written once and stay the same Julia object (`===`) after loading.

```@docs; canonical=false
save_setup
load_setup
```

## Example

```julia
using BeamletOptics

m1 = SquarePlanoMirror2D(0.0254)
m2 = SquarePlanoMirror2D(0.0254)
translate3d!(m2, [0, 0.1, 0])
system = System([m1, m2])
beam = GaussianBeamlet([0, -0.1, 0], [0, 1.0, 0], 635e-9, 1e-4)

save_setup("setup.bmo", system => beam)

setup = load_setup("setup.bmo")
loaded_system, loaded_beam = only(setup.pairs)
```

The archive is a zip file with a `setup.toml` description and an `assets` folder for
binary data (e.g. mesh vertices and faces); identical assets are stored only once.

## Compatibility

A setup file can only be loaded by a BeamletOptics release from the same breaking series
that wrote it (same major version, or same minor version for `0.x` releases) and that is
not older than the file itself; see [`load_setup`](@ref) above. Any other combination
raises an error naming the BeamletOptics version needed to open the file.

Only the untraced state of a setup is stored: sources are written as they were constructed,
not as traced, so a [`Detector`](@ref)'s `hits` are never part of the file, and a loaded
source has to be traced again with [`solve_system!`](@ref).

## Refractive indices

An object's refractive index (see [`RefractiveIndex`](@ref)) is stored differently
depending on its type:

- a plain number becomes a [`ConstantRefractiveIndex`](@ref);
- [`SellmeierEquation`](@ref) and [`DiscreteRefractiveIndex`](@ref) are stored with their
  parameters;
- any other callable, e.g. a closure or a named function, must first be registered with
  [`register_material!`](@ref) under a name. Only that name is written to the file, so the
  same name must be registered again — with the same refractive index — before the file is
  loaded.

```@docs; canonical=false
ConstantRefractiveIndex
register_material!
```

## Extending storage to custom types

A custom [`BeamletOptics.AbstractObject`](@ref) subtype or beam source becomes storable by
implementing `BeamletOptics.to_storage`/`BeamletOptics.from_storage` and registering a tag
with `BeamletOptics.register_storage_type!`. Packages that define storable types should call
`register_storage_type!` from their `__init__` function, so the registration happens once
per Julia session regardless of load order.

Binary data that does not belong in `setup.toml` (e.g. mesh vertices and faces) is written
to the archive's `assets` folder with `BeamletOptics.write_asset!` and read back with
`BeamletOptics.read_asset`; `BeamletOptics.encode_array`/`BeamletOptics.decode_array` do
this automatically for arrays of `Float64`, `Int64` or `ComplexF64` above a size threshold,
keeping small arrays inline in the TOML instead.

```@docs; canonical=false
BeamletOptics.register_storage_type!
BeamletOptics.to_storage
BeamletOptics.from_storage
BeamletOptics.write_asset!
BeamletOptics.read_asset
BeamletOptics.encode_array
```
