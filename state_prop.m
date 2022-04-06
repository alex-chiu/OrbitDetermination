function [dY] = state_prop(time, Y, t0_TAI, TwoB, J2, D, AMat)
% full_prop : Full propagation function with perturbations
%
% INPUTS
%
% time ------ Time value for propagation, in s
%
% Y --------- 132-by-1 state vector from previous time step
%
% t0_TAI ---- Initial Julian date atomic time, in days
%
% C --------- Structure containing all relevant constants
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
% Last Edited: 3/29/2022
%+==============================================================================+

dY = zeros(132, 1);

% Extract values
R_ECI = Y(1:3);
V_ECI = Y(4:6);
Phi = reshape(Y(12:end), 11, 11);

% Time calculations
t_TAI = t0_TAI + time / 86400; 
t_TT = t_TAI + 32.184 / 86400;

% Find acceleration
A_ECI = TwoB(R_ECI(1), R_ECI(2), R_ECI(3)) + ...
        J2(R_ECI(1), R_ECI(2), R_ECI(3)) + ...
        D(Y(7), V_ECI(1), V_ECI(2), V_ECI(3),   ...
                 R_ECI(1), R_ECI(2), R_ECI(3));

% Populate dY
dY(1:3) = V_ECI;
dY(4:6) = A_ECI;

A = AMat(Y(7), V_ECI(1), V_ECI(2), V_ECI(3), R_ECI(1), R_ECI(2), R_ECI(3));
dPhi = A * Phi;
dY(12:end) = reshape(dPhi, 121, 1);