# OrbitDetermination
Code for Methods of Orbit Determination

Main Program - Final_Project_Main.m

get_symb_func.m
Calculates symbolic functions for finding A and \tilde{H} Matrices, and G Function. Returns a set of function handles.

Tips (4/7/2022)

Assumptions:
- ECEF to ECI: 10 cm
- Force Model: 20 x 20 EGM96, SRP, Drag, 3rd Body
- Estimated Parameters: Position, Velocity in ECI
- Process Noise: 1 x 10^{-12} km/s^2
- A Priori State
- A Priori Covariance: Position: 1 km, Velocity: 10 m/s
- Estimation Strategy: EKF
- Q [1 x 10^{-12}, 1 x 10^{-9} 1 x 10^{-12}]

Model light-time correction
- Propagated shift
- Ground station aberration

Have flags for various cases:
- Case of interest
- SRP (Cannonball, Facet, Attitude)
- Drag (Area - x, Facet)

Test Certain Type of Measurement - Turn up measurement noise for the measurement you don't want information from.

TODO:
- Work on Dynamics
    - 20 x 20 gravity model
    - Refine drag
    - Add SRP 
    - Add 3rd body
- Check filter and filter equations
- Lighttime correction