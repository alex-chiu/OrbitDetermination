%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASE 389.4 Methods of Orbit Determination
% Final Project
% Author: Alex Chiu (ac68767)
% Last Edited: 04/05/2022
% Summary: Main Script for Final Project
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear Workspace
clc; clear all; close all; format long g;

% Start Timer
tic;

% Load Constants and Initial Conditions
constants_script;
initial_conditions;

% Get Function Handles for Various Matrices
[getAMat, getHTildeMat, getGVec] = get_symb_func(C);

A_0 = getAMat(2.0, ...                                         % C_d
              IC.V_ECI_0(1), IC.V_ECI_0(2), IC.V_ECI_0(3), ... % Velocity
              IC.R_ECI_0(1), IC.R_ECI_0(2), IC.R_ECI_0(3));    % Position

A_sol = load("..\HW 5\Data\A_t0.mat").A;
A_diff = abs((A_0(1:7, 1:7) - A_sol) ./ A_sol);
disp(A_diff);

% End Timer
toc;