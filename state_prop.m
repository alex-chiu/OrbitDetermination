function [dY] = state_prop(time, Y, t0_UTC, getAMat, getAcc)
% full_prop : Full propagation function with perturbations
%
% INPUTS
%
% time ------ Time value for propagation, in s
%
% Y --------- 132-by-1 state vector from previous time step
%
% t0_UTC ---- Initial Julian date in UTC
%
% getAMat --- Function handle for A matrix
%
% getAcc ---- Function handle for acceleration calculations
%
% OUTPUTS
%
% dY -------- 132-by-1 time derivative vector
% 
%+------------------------------------------------------------------------------+
% References: Lecture Notes 
%
% Author: Alex Chiu
%
% Last Edited: 4/8/2022
%+==============================================================================+

dY = zeros(132, 1);

% Extract values
R_ECI = Y(1:3);
V_ECI = Y(4:6);
Phi = reshape(Y(12:end), 11, 11);

% Time calculations
t0_TAI = t0_UTC + 37 / 86400;
t_TT = t0_TAI + (time + 32.184) / 86400;

% 20x20 EGM96 Model
% R_ECEF = QMat' * R_ECI;
% [G20x, G20y, G20z] = gravitysphericalharmonic(R_ECEF', 'EGM96', 20, 'None');
% AccG20 = [G20x; G20y; G20z] / 1000;

% Find acceleration
A_ECI = getAcc(Y(7), V_ECI(1), V_ECI(2), V_ECI(3), R_ECI(1), R_ECI(2), R_ECI(3));
% A_ECI = getAcc(Y(7), V_ECI(1), V_ECI(2), V_ECI(3), R_ECI(1), R_ECI(2), R_ECI(3)) + AccG20;

% Populate dY
dY(1:3) = V_ECI;
dY(4:6) = A_ECI;

A = getAMat(Y(7), V_ECI(1), V_ECI(2), V_ECI(3), R_ECI(1), R_ECI(2), R_ECI(3));
dPhi = A * Phi;
dY(12:end) = reshape(dPhi, 121, 1);