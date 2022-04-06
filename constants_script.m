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
C.Area = 6;

% Atmospheric Parameters
% Reference height for atmosphere model, in km
C.H = 88.667; 
% Base density of air, in kg/m^3
C.rho_0 = 3.614e-13;
% Base radius, in km;
C.r_0 = 700 + C.R_e; 