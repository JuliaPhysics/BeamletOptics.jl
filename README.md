# SCDI Sim

This project implements a forward model to simulate different aspects of a single crystal dispersion interferometer (SCDI) setup.

## Features 

- Complex ray tracing
    - Utilities
        - [x] basic translation/rotation tools
        - [x] basic STL data input
        - [ ] STL vertices/face compression?
    - Rays        
        - [x] basic intersection testing
        - [x] basic propagation routine over n elements
        - [ ] basic optical elements (lens/mirror)
            - [x] reflection
            - [ ] refraction
- Optics
    - [ ] phase front via Gaussian beams
    - [ ] phase shift due to ref. index change
    - [ ] optical intensity distribution on detector
- Mechanics
    - [ ] vibration of optical elements
- Plotting
    - [ ] automatic plot updates using Makie Observable/Buffer
- Test coverage
    - [ ] translation and x,y,z-rotation tests
    - [ ] Möller-Trumbore-algorithm test
    - [ ] generic mesh intersection test

## Known bugs

- [x] wrong intersection for some geometries/angles?
    * fixed typo in Möller-Trumbore-algorithm
- [ ] align3d rotation incorrect if start = -target
- [ ] red warntype for trace_system
    * even though SCDI.intersect3d(object::Geometry, ray::Ray) is type stable?
- [x] normal calculated from face in orthogonal3d(object::Geometry, fID::Int) points in the wrong dir?

## Comments

* the struct that contains the ray data (pos,dir,len) is currently mutable in order to change the length
    * ray length is initialized with `Inf`
    * ray length will be changed to `t` if an intersection occurs
    * thus requires mutable struct
    * is there a better (syntactial) way to do this, i.e. with immutable?
* "Inheritance" is solved via composition
    * macro injects necessary fields into struct
    * macro ensures all wrapper functions are defined
