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
        - [ ] basic propagation routine over n elements
        - [ ] basic optical elements (lens/mirror)
            - [ ] reflection
            - [ ] refraction
- Optics
    - [ ] phase front via Gaussian beams
    - [ ] phase shift due to ref. index change
    - [ ] optical intensity distribution on detector
- Mechanics
    - [ ] vibration of optical elements
- Test coverage
    - [ ] translation and x,y,z-rotation tests
    - [ ] Möller-Trumbore-algorithm test
    - [ ] generic mesh intersection test

## Known bugs

- [x] wrong intersection for some geometries/angles?
    * fixed typo in Möller-Trumbore-algorithm
- [ ] align3d rotation incorrect if start = -target

## Comments

* the struct that contains the ray data (pos,dir,len) is currently mutable in order to change the length
    * ray length is initialized with `Inf`
    * ray length will be changed to `t` if an intersection occurs
    * thus requires mutable struct
    * is there a better (syntactial) way to do this, i.e. with immutable?