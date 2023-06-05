# SCDI Sim

This project implements a forward model to simulate different aspects of a single crystal dispersion interferometer (SCDI) setup.

## Documentation

Documentation can be found here: https://optical-air-data.pages.gitlab.dlr.de/dispersionsinterferometer/scdi-sim/

## Features 

- Complex ray tracing
    - Utilities
        - [x] basic translation/rotation tools
        - [x] basic STL data input
        - [x] reset functionality for rotation/translation
    - Rays        
        - [x] basic intersection testing
        - [x] basic propagation routine over n elements
        - [x] basic optical elements (lens/mirror)
            - [x] reflection
            - [x] refraction
        - [x] Gaussian beam struct
            - [ ] implement better retracing based on chief ray
            - [ ] simple astigmatism via 5, 7, 9-ray approach
- Optics
    - [x] Lens types using mathematical surfaces
    - [ ] Photodetector
        - [ ] optical intensity distribution on detector
        - [ ] model interference of Gaussian beams
        - [ ] phase front via Gaussian beams
    - [ ] phase shift due to ref. index change
- Mechanics
    - [ ] vibration of optical elements
- Plotting
    - [ ] automatic plot updates using Makie Observable/Buffer
- Test coverage
    - [x] continuous integration pipeline by O. Kliebisch
    - [x] translation and x,y,z-rotation tests
    - [x] Möller-Trumbore-algorithm test
    - [x] generic mesh intersection test
    - [x] solve_system()/trace_system() test
    - [ ] interact() tests 
    - [ ] more Lens/Meshlikes tests

## Known bugs

- [ ] align3d rotation incorrect if start = -target
- [ ] @code_warntype for interact(Lens, Beam)
- [ ] change names for functions that return matrices
- [x] weird results for intersect(::Sphere, ::Ray)
    - fixed by O.K.
- [ ] weird results for high level-of-detail meshes

## Comments

- none
