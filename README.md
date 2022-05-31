# SCDI Sim

This project implements a forward model to simulate different aspects of a single crystal dispersion interferometer (SCDI) setup.

## Features 

- Complex ray tracing
    - Utilities
        - [] basic translation/rotation tools
        - [] basic STL data input
        - [] STL vertices/face compression?
    - Rays        
        - [] basic intersection testing
        - [] basic optical elements (lens/mirror)
        - [] reflection
        - [] refraction
- Optics
    - [] phase front via Gaussian beams
    - [] phase shift due to ref. index change
    - [] optical intensity on detector
- Mechanics
    - [] vibration of optical elements
- Test coverage
    - [] translation and x,y,z-rotation tests
    - [] Möller-Trumbore-algorithm test
    - [] generic mesh intersection test

## Known bugs

- [] align3d rotation incorrect if $v_{\mathrm{start}} = -v_{\mathrm{target}}$ 