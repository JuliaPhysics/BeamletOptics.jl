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
            - [x] implement better retracing based on chief ray
            - [ ] simple astigmatism via 5, 7, 9-ray approach
- Optics
    - [x] Lens types using mathematical surfaces
    - [x] Interference
        - [x] optical intensity distribution on detector
        - [x] model interference of Gaussian beams
        - [x] phase front via Gaussian beams
        - [ ] model coherence length via loss of contrast
    - [ ] phase shift due to ref. index change
    - [ ] Multi-body container type, i.e. telescope, for easy kinematics
- Mechanics
    - [x] element group handling
    - [ ] vibration of optical elements
- Plotting
    - [ ] automatic plot updates using Makie Observable/Buffer
- Test coverage
    - [x] continuous integration pipeline by O. Kliebisch
    - [x] translation and x,y,z-rotation tests
    - [x] Möller-Trumbore-algorithm test
    - [x] generic mesh intersection test
    - [x] solve_system()/trace_system() test
    - [x] interact() tests 
    - [x] more Lens/Meshlikes tests

## Known bugs

- [ ] interact3d for AbstractRefractiveOptic/Abstract Ray: n2 not correctly set in case of TIR
- [ ] inconsistent use of nm, m, mm (i.e. in Ray constructor)
- [x] align3d rotation incorrect if start = -target
- [ ] @code_warntype for interact(Lens, Beam)
- [ ] change names for functions that return matrices
- [ ] weird results for high level-of-detail meshes
- [x] SDFs and meshes rotate in opposite directions
- [x] GaussianBeamlet efield calculation fails at beam waist
- [ ] point_on_beam calculation can be incorrect at optical surfaces
- [ ] intersection calculation can be incorrect if two objects "touch"
- [ ] inconsistent rotation behavior for shapes of type Mesh 

## Comments

- none
