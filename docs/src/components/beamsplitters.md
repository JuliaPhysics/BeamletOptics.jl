# Beamsplitters

Beamsplitters are used to split beams of light, enabling the separation of an incoming beam into reflected and transmitted parts. In this package, beamsplitters are implemented via the [`BeamletOptics.AbstractBeamsplitter`](@ref) interface. This type loosely defines the interaction logic used for the tracing and retracing of optical systems that incorporate these devices. When a beam is
split the new beams are referred to as the `children` of the `parent` beam. This package uses the `AbstractTrees.jl` interface in order to represent the splitting of [`BeamletOptics.AbstractBeam`](@ref)s.
For concrete implementations of splitters, skip to the [Rectangular plate beamsplitter](@ref) or [Cube beamsplitter](@ref) section.

!!! info "Child order"
    In order to ensure consistency the following definition for the order of appended child beams is used: The **transmitted beam is appended first** and the reflected beam is appended second.
    The beam path of the `parent` beam stops at the splitting interface, i.e. a `nothing` optical interaction occurs.


## Thin beamsplitter

This model acts as a quasi-coating, representing a very thin (zero-thickness) layer that directly splits the beamlets. This type is used for testing purposes mainly and to build composite objects. For optical setups it is recommended to use one of the splitter types listed below. One of two constructors can be used in order to spawn thin splitters:

- [`ThinBeamsplitter`](@ref)
- [`RoundThinBeamsplitter`](@ref)

!!! info "Beamsplitter phase jump"
    If you want to learn more about how beam splitting is modeled, especially in the context of interferometer simulations, refer to the [`ThinBeamsplitter`](@ref) type docs and the specific [`BeamletOptics.interact3d`](@ref) documentation.

## Plate beamsplitters

Plate beamsplitters provide a more sophisticated model by incorporating a substrate with finite thickness and [`BeamletOptics.RefractiveIndex`](@ref). The splitter coating is modeled via a [`ThinBeamsplitter`](@ref) placed flush onto a single face of the substrate. This allows for more realistic simulation of refractive effects such as e.g. beam path displacement and optical path differences caused by the substrate's geometry. Below two concrete implementations are showcased.

### Rectangular plate beamsplitter

The [`RectangularPlateBeamsplitter`](@ref) represents a planar, rectangular substrate with a partially reflective coating. 

```@docs; canonical=false
RectangularPlateBeamsplitter(::Real, ::Real, ::Real, ::BeamletOptics.RefractiveIndex)
```

A common application involves that this type of beamsplitter is paired with a compensator plate. Below an exemplary illustration of such a setup is shown, where the beamsplitter reflects part of the incoming beam perpendiculary. A [`RectangularCompensatorPlate`](@ref) ensures that the parallel path offset is corrected.

```@eval
file_dir = joinpath(@__DIR__, "..", "assets", "bs_assets")

Base.include(@__MODULE__, joinpath(file_dir, "pbs_showcase.jl"))

take_screenshot("pbs_showcase.png", system, beam; size=(600, 400), view=pbs_view, px_per_unit=8, optional_rendering_fct=test)
```

![Plate beamsplitter showcase](pbs_showcase.png)

The depicted system consist of a rectangular beamsplitter (e.g. [BSW26R](https://www.thorlabs.com/thorproduct.cfm?partnumber=BSW26R)) in combination with a compensator plate (e.g. [BCP42R](https://www.thorlabs.com/thorproduct.cfm?partnumber=BCP42R)). They are both mounted in [KM2536](https://www.thorlabs.com/thorproduct.cfm?partnumber=KM2536) kinematic mounts. The splitter substrate thickness is exaggerated for the purpose of illustration.

### Round plate beamsplitter

This variant uses a circular substrate. While components of this type often feature a wedge angle to avoid "ghost beams" in practice, e.g. the [BSW26](https://www.thorlabs.com/thorproduct.cfm?partnumber=BSW26), this is not modeled here.

```@docs; canonical=false
RoundPlateBeamsplitter(::Real, ::Real, ::BeamletOptics.RefractiveIndex)
```

## Cube beamsplitter

The [`CubeBeamsplitter`](@ref) is composed of two [`RightAnglePrism`](@ref)s with a partially reflective interface at their internal joint. The splitting interface is represented by a [`ThinBeamsplitter`](@ref). As with other beamsplitters in this package, the cube beamsplitter uses scalar reflection and transmission coefficients.

```@docs; canonical=false
CubeBeamsplitter(::Real, ::BeamletOptics.RefractiveIndex)
```

```@eval
file_dir = joinpath(@__DIR__, "..", "assets", "bs_assets")

Base.include(@__MODULE__, joinpath(file_dir, "cbs_showcase.jl"))

take_screenshot("cbs_showcase.png", system, beam; size=(600, 400), view=mzi_view, px_per_unit=8)
```

A classic application of cube beamsplitters is in Mach–Zehnder interferometers, where two beamsplitters combine with additional mirrors to form two optical paths that later recombine. The figure below shows a rudimentary Mach–Zehnder arrangement using two cube beamsplitters and two [`RightAnglePrismMirror`](@ref)s:

![Cube beamsplitter showcase](cbs_showcase.png)

The system is made up of two cube beamsplitters (e.g. [BS013](https://www.thorlabs.com/thorproduct.cfm?partnumber=BS013)) and two mirrors mounted in [KM100PM/M](https://www.thorlabs.com/thorproduct.cfm?partnumber=KM100PM/M) prism mounts with [PM4](https://www.thorlabs.com/thorproduct.cfm?partnumber=PM4) clamping arms. The [`GaussianBeamlet`](@ref) enters from the bottom left and exits at the top right. Note that due to the splitting, two beams exit the systems. Depending on the relative phase, they mutually interfere. Total optical power is conserved.

!!! tip "Interferometer tutorial"
    Refer to the [Michelson interferometer](@ref) tutorial for a detailed showcase featuring beam splitters and other components.