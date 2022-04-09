% Loads all necessary constants for orbit determination process

% Earth Parameters
% Gravitational constant of Earth, in km^3/s^2
C.mu_e = 398600.4415;
% Radius of Earth, in km
C.R_e = 6378.1363;
% Earth J2 Constant
C.J2 = 0.00108248;
% Rotation rate of Earth, in rad/s
C.w_e = 7.292115146706979e-5;

% Spacecraft Parameters
% Spacecraft mass, in kg;
C.m = 2000;

% Spacecraft areas, in m^2;
C.area = 6;

% Atmospheric Parameters
% Reference height for atmosphere model, in km
C.H = 88.667; 
% Base density of air, in kg/m^3
C.rho_0 = 3.614e-13;
% Base radius, in km;
C.r_0 = 700 + C.R_e; 
% Drag Coefficient
C.C_d = 2.0;

% Station Parameters
% Positions in ECEF (ITRF), in km
C.station_pos = [-6143.584, 1907.295, 2390.310;
                 1364.250, 6030.810, -5564.341;
                 1033.743, -817.119, 1994.578];

% ODE Parameters
C.options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);