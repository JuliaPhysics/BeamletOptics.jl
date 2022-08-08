# SCDI Sim

This project implements a forward model to simulate different aspects of a single crystal dispersion interferometer (SCDI) setup.

## Features 

- Complex ray tracing
    - Utilities
        - [x] basic translation/rotation tools
        - [x] basic STL data input
        - [ ] STL vertices/face compression?
        - [x] reset functionaltity for rotation/translation
    - Rays        
        - [x] basic intersection testing
        - [x] basic propagation routine over n elements
        - [x] basic optical elements (lens/mirror)
            - [x] reflection
            - [x] refraction
            - [ ] documentation
        - [x] wavelength field for Ray
        - [ ] Gaussian beam struct
- Optics
    - [ ] phase front via Gaussian beams
    - [ ] phase shift due to ref. index change
    - [ ] optical intensity distribution on detector
    - [ ] simple astigmatism via 5, 7, 9-ray approach
- Mechanics
    - [ ] vibration of optical elements
- Plotting
    - [ ] automatic plot updates using Makie Observable/Buffer
- Test coverage
    - [x] translation and x,y,z-rotation tests
    - [ ] Möller-Trumbore-algorithm test
    - [ ] generic mesh intersection test

## Known bugs

- [x] wrong intersection for some geometries/angles?
    * fixed typo in Möller-Trumbore-algorithm
- [ ] align3d rotation incorrect if start = -target
- [x] @code_warntype for trace_system
    * fixed by O. Kliebisch by stabilizing intersect3d
- [x] normal calculated from face in orthogonal3d(object::Geometry, fID::Int) points in the wrong dir?
    * fixed typo in reflection calculation
- [ ] @code_warntype for interact(Lens, Beam)
- [ ] change names for functions that return matrices

## Comments

* the struct that contains the ray data (pos,dir,len) is currently mutable in order to change the length
    * ray length is initialized with `Inf`
    * ray length will be changed to `t` if an intersection occurs
    * thus requires mutable struct
    * is there a better (syntactial) way to do this, i.e. with immutable?
* Möller-Trumbore-algorithm
    * det < kepsilon condition: use of abs(det) or not?
    * original paper does not feature abs()
    * modern implementations do use abs(det)
    * errors during testing if used without abs()!
* type structure proposal:
    * AbstractEntity
        * AbstractNoMeshPlaceholder
        * AbstractMesh
            * AbstractBoundedMesh
                * Prism
                * Mirror
* replace Matrix representation of faces and vertices by Face type?
    * only complicates plotting etc. and not really worth it
