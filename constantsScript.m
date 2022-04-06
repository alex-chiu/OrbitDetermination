% A script that loads constants used for OD

% Gravitational constant of Earth, in km^3/s^2
C.mu_e = 398600.4415;
% Radius of Earth, in km
C.R_e = 6378.1363;

% J2 Constant
C.J2 = 0.00108248;

% Spacecraft mass, in kg;
C.m = 2000;
% Area of spacecraft, in m^2
C.A_sp = 6.0;
% Rotation rate of Earth, in rad/s
C.w_e = 7.292115146706979e-5;
% Reference height for atmosphere model, in km
C.H = 88.667; 
% Base density of air, in kg/m^3
C.rho0 = 3.614e-13;
% Base radius, in km;
C.r0 = 700 + C.R_e; 
% Coefficient of Drag
C.C_d = 100 / 94 * 1.88;
% Solar Coefficient
C.C_solar = 0;

% Measurement Errors
C.sigmaR = 0.005; % Range (km)
C.sigmaRR = 0.000001; % Range (km/s)

C.start_time = 2458150.70833; % JD

% Station Positions
C.station_pos = [-6143.584, 1907.295, 2390.310;
                 1364.250, 6030.810, -5564.341;
                 1033.743, -817.119, 1994.578];