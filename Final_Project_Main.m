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

% Flag for testing initial conditions
test_init = false;

if (test_init)
    [PMat, NMat, SMat, WMat] = get_ECEF_ECI(2458150.70833, 0.1965014, 0.015361, 0.288259); %#ok<*UNRCH> 
    PNSMat = PMat * NMat * SMat;
    R_S_ECEF = C.station_pos(:, 1);
    
    A_0 = getAMat(2.0, ...                                         % C_d
                  IC.V_ECI_0(1), IC.V_ECI_0(2), IC.V_ECI_0(3), ... % Velocity
                  IC.R_ECI_0(1), IC.R_ECI_0(2), IC.R_ECI_0(3));    % Position
    
    A_sol = load("..\HW 5\Data\A_t0.mat").A;
    A_diff = abs((A_0(1:7, 1:7) - A_sol) ./ A_sol);
    disp(A_diff);
    
    H_0 = getHTildeMat(PNSMat(1, 1), PNSMat(1, 2), PNSMat(1, 3), ...
                       PNSMat(2, 1), PNSMat(2, 2), PNSMat(2, 3), ...
                       PNSMat(3, 1), PNSMat(3, 2), PNSMat(3, 3), ...
                       R_S_ECEF(1), R_S_ECEF(2), R_S_ECEF(3), ...       % Station Position in ECEF
                       WMat(1, 1), WMat(1, 2), WMat(1, 3), ...
                       WMat(2, 1), WMat(2, 2), WMat(2, 3), ...
                       WMat(3, 1), WMat(3, 2), WMat(3, 3), ...
                       IC.V_ECI_0(1), IC.V_ECI_0(2), IC.V_ECI_0(3), ... % Velocity
                       IC.R_ECI_0(1), IC.R_ECI_0(2), IC.R_ECI_0(3));    % Position
    
    H_sol = load("..\HW 5\Data\H_Tilde_t0.mat").H_TILDA;
    H_diff = abs((H_0(:, 1:7) - H_sol) ./ H_sol);
    disp(H_diff);
end

% Load Measurements
MeasMat = load("..\HW 5\Data\LEO_DATA_Apparent.mat").LEO_DATA_Apparent;

% Run Extended Kalman Filter


% End Timer
toc;