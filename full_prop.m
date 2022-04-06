function [dY] = full_prop(time, Y, t0_TAI, TwoB, J2, D, AMat, C)
% full_prop : Full propagation function with perturbations
%
% INPUTS
%
% time ------ Time value for propagation, in s
%
% Y --------- 55-by-1 state vector from previous time step
%
% t0_TAI ---- Initial Julian date atomic time, in days
%
% C --------- Structure containing all relevant constants
%
% OUTPUTS
%
% a --------- 3-by-1 acceleration vector
% 
%+------------------------------------------------------------------------------+
% References: Lecture Notes 
%
% Author: Alex Chiu
%
% Last Edited: 3/29/2022
%+==============================================================================+

dY = zeros(55, 1);

% Extract values
R_ECI = Y(1:3);
V_ECI = Y(4:6);
Phi = reshape(Y(7:end), 7, 7);

% Time calculations
t_TAI = t0_TAI + time / 86400; 
t_TT = t_TAI + 32.184 / 86400;

% Find gravity acceleration

A_ECI = TwoB(R_ECI(1), R_ECI(2), R_ECI(3)) + ...
        J2(R_ECI(1), R_ECI(2), R_ECI(3)) + ...
        D(C.C_d, V_ECI(1), V_ECI(2), V_ECI(3),   ...
                 R_ECI(1), R_ECI(2), R_ECI(3));

% Populate dY
dY(1:3) = V_ECI;
dY(4:6) = A_ECI;

A = AMat(C.C_d, V_ECI(1), V_ECI(2), V_ECI(3), R_ECI(1), R_ECI(2), R_ECI(3));

dPhi = A(1:7, 1:7) * Phi;
dY(7:end) = reshape(dPhi, 49, 1);