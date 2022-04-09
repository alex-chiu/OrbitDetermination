function [getAMat, getHTildeMat, getGVec, getAcc] = get_symb_func(C, State_Flag)
% full_prop : Full propagation function with perturbations
%
% INPUTS
%
% C -------------- Structure containing all necessary constants
%
% State_Flag ----- Flag for number of states to use
%
% OUTPUTS
%
% getAMat -------- Function handle for calculating the A matrix. 
%
% getHTildeMat --- Function handle for calculating the HTilde matrix.
%
% getGVec -------- Function handle for calculating the G vector.
% 
%+------------------------------------------------------------------------------+
% References: Lecture Notes 
%
% Author: Alex Chiu
%
% Last Edited: 4/5/2022
%+==============================================================================+

% Define Storage Matrices
A = sym("A", [11 11]);
HTilde = sym("HTilde", [2 11]);
G = sym("G", [2 1]);

% Define State Variables
syms x y z vx vy vz C_d C_solar b1 b2 b3;

if (State_Flag == 6)
    C_d = C.C_d;
end

% Populate Helper Vectors
Pos = [x; y; z];
Vel = [vx; vy; vz];

% Find Two Body Acceleration
Acc_2B = -C.mu_e / norm(Pos)^3 * Pos;

% Find J2 Acceleration
Acc_J2 = sym('Acc_J2', [3 1]);
Acc_J2_Mag = 3 * C.mu_e * C.J2 * C.R_e^2 / (2 * norm(Pos)^5);
Acc_J2(1) = Acc_J2_Mag * x * (5 * (z/norm(Pos))^2 - 1);
Acc_J2(2) = Acc_J2_Mag * y * (5 * (z/norm(Pos))^2 - 1);
Acc_J2(3) = Acc_J2_Mag * z * (5 * (z/norm(Pos))^2 - 3);

% Find Drag Acceleration
V_Rel = [vx + C.w_e*y; vy - C.w_e*x; vz];
Acc_Drag = -1/2 * C_d * C.area / C.m * C.rho_0 * ...
           exp(-(norm(Pos) - C.r_0) / C.H) * 1000 * norm(V_Rel) * V_Rel;

% Find SRP Acceleration
% TODO

% Combine Accelerations
Acc = Acc_2B + Acc_J2 + Acc_Drag;

% Find Range and Range Rate
R_S_ECEF = sym('R_S_ECEF', [3 1]);
WMat = sym('WMat', [3 3]);
PNSMat = sym('PNSMat', [3 3]);
R_S_PEF = WMat * R_S_ECEF;
R_S_ECI = PNSMat * R_S_PEF;
V_S_ECI = PNSMat * cross([0; 0; C.w_e], R_S_PEF);

Range = norm(Pos - R_S_ECI);
Range_Rate = dot(Pos - R_S_ECI, Vel - V_S_ECI) / Range;

% Form Full Vectors
X = [Pos; Vel; C_d; C_solar; b1; b2; b3];
F = [Vel; Acc; 0; 0; 0; 0; 0];
Z = [Range; Range_Rate];

% Find Partial Derivatives for Matrices
for i = 1:11
    for j = 1:11
        A(i, j) = diff(F(i), X(j));
    end
    for k = 1:2
        HTilde(k, i) = diff(Z(k), X(i));
    end
end

% Convert Matrices to Function Handles
getAMat = matlabFunction(A);
getHTildeMat = matlabFunction(HTilde);
getGVec = matlabFunction(Z);
getAcc = matlabFunction(Acc);
end