function [getAMat, getHTildeMat, getGVec] = get_symb_func(C)
% full_prop : Full propagation function with perturbations
%
% INPUTS
%
% C -------------- Structure containing all necessary constants
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
Acc_Drag = -1/2 * C_d * C.Area / C.m * C.rho_0 * ...
           exp(-(norm(Pos) - C.r_0) / C.H) * 1000 * norm(V_Rel) * V_Rel;

% Combine Accelerations
Acc = Acc_2B + Acc_J2 + Acc_Drag;

% Form Full Vectors
X = [Pos; Vel; C_d; C_solar; b1; b2; b3];
F = [Vel; Acc; 0; 0; 0; 0; 0];

% Find Partial Derivatives for A Matrix
for i = 1:11
    for j = 1:11
        A(i, j) = diff(F(i), X(j));
    end
end

% Convert Matrices to Function Handles
getAMat = matlabFunction(A);
getHTildeMat = matlabFunction(HTilde);
getGVec = matlabFunction(G);
end