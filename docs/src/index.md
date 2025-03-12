# BeamletOptics

Building optical setups in a laboratory environment -- for instance a laser interferometer -- is a common task for optical engineers and physicists. This package is intended to provide a simulation environment in which the user can quickly analyze and layout simple optical components like lenses or beamsplitters before committing to a breadboard setup.  

## What is the purpose of this package

This package mainly tries to provide a simple Gaussian beamlet propagation tool for coherent, monochromatic and directed light sources. It also offers a convenient kinematic API that allows for the easy placing of optical elements and straight-forward simulation of moving or vibrating components. 

For this purpose, the package implements a traditional ray tracing solver. This forms the backbone of the Gaussian beamlet tracing scheme that has been implemented to model the propagation of laser beams.

!!! info "What this package is not"
    This package does not include tools for optimizing optical systems, such as fine-tuning lens surfaces to minimize specific aberrations in multi-lens setups. Instead, the package is designed as a digital laboratory where you can play around with stuff before buying it.

## Features list

- Hybrid sequential and non-sequential 3D ray tracing without paraxial approximation
- TEM₀₀ [Gaussian beamlet](@ref) models
- Various optical components
    - [Mirrors](@ref)
    - [Lenses](@ref)
    - [Beamsplitters](@ref)
    - [Detectors](@ref)
- Surface-like modeling of rotationally symmetrical lens systems
- Extendable [API design](@ref) for the implementation of custom optical interactions
- Easy visualization via the Makie package

## Getting started

Refer to the sections below if this is your first time using this package.

### Installation

!!! warning
    This package requires Julia ≥ 1.9.4

You can add this package to your project by entering the package manager (press `]` in the REPL) and typing `add BeamletOptics`. It is also recommended that you `add GLMakie`. You can include this package into your current scope via `using BeamletOptics`. If a Makie version is loaded before or after the inclusion of this package, the extension provided as part of this package will enable additional visualization functions. 

### Tutorials

Follow the [Beam expander](@ref) and [Michelson interferometer](@ref) tutorials for a quick start into the package interface.

## Citation and license

The BeamletOptics package is made available under the MIT license. If you use this package for your research, we encourage you to cite it. For your convenience, a BibTeX entry is provided in the **FIXME** file.

## Similar packages

A variety of packages and tools exist that implement similar approaches or offer optics modeling capabilities. Within the Julia ecosystem, the following packages need to be mentioned:

- [OpticSim.jl](https://github.com/brianguenter/OpticSim.jl)
- [ABCDMatrixOptics.jl](https://github.com/JuliaPhysics/ABCDMatrixOptics.jl)
- [WaveOpticsPropagation.jl](https://github.com/JuliaPhysics/WaveOpticsPropagation.jl)

More broadly speaking, have a look at these packages as well:

- [DynamicalBilliards.jl](https://github.com/JuliaDynamics/DynamicalBilliards.jl)
- [RayTracer.jl](https://github.com/avik-pal/RayTracer.jl)

There also exists a plethora of commercial and non-commercial simulation frameworks outside of the Julia ecosystem. For specific examples regarding the beamlet method used in this package, refer to the [Complex ray tracing](@ref) section. 

## Development roadmap

In order to warrant a 1.0.0 release tag, the following features will need to be implemented. This definition is arbitrary. The exact timeframe for this development effort is not specified, but will be on the order of 1-2 years.

- 🔳 TODO
- 🟩 WIP
- ✅ DONE

### Additional features

- 🔳 Beamlet tracing
    - 🟩 Implementation of the full polarized astigmatic Gaussian beamlet formalism
        - refer to Worku et al. (2017/2020), Greynolds (1985)
    - 🔳 Support for (trivial) forms of 2D-field decomposition, e.g. tophat
    - 🔳 Modeling of the coherence length and contrast influence
- 🔳 Components
    - 🟩 Polarizing optics (based on Jones formalism)
        - 🔳 Wave plates (λ/2 and λ/4)
        - 🔳 Retarders
        - 🔳 Faraday rotators
    - 🔳 Isolators
    - 🔳 Modulators
        - 🔳 AOM
        - 🔳 EOM
- 🔳 Visualization
    - 🔳 automatic Makie plot updates (e.g. some form of "interactive" mode)
    - 🔳 Better Lens surface plots based on multiple dispatch