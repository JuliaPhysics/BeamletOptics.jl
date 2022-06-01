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
        - [ ] basic optical elements (lens/mirror)
        - [ ] reflection
        - [ ] refraction
- Optics
    - [ ] phase front via Gaussian beams
    - [ ] phase shift due to ref. index change
    - [ ] optical intensity on detector
- Mechanics
    - [ ] vibration of optical elements
- Test coverage
    - [ ] translation and x,y,z-rotation tests
    - [ ] Möller-Trumbore-algorithm test
    - [ ] generic mesh intersection test

## Known bugs

- [x] wrong intersection for some geometries/angles?
- [ ] align3d rotation incorrect if start = -target 