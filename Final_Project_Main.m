%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASE 389.4 Methods of Orbit Determination
% Final Project
% Author: Alex Chiu (ac68767)
% Last Edited: 04/09/2022
% Summary: Main Script for Final Project
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO:
% - Upgrade to 20 x 20 EGM96 Gravity
% - Add SRP
% - Upgrade Drag Model
% - Add 3rd Body Effects
% - Add Light-Time Effects

% Clear Workspace
clc; clear all; close all; format long g;

% Start Timer
tic;

% Load Constants and Initial Conditions
constants_script; initial_conditions;

% Get Function Handles for Various Matrices
[getAMat, getHTildeMat, getGVec, getAcc] = get_symb_func(C, 6);

% Load Measurements
MeasMat = load(".\Data\LEO_DATA_Apparent-1.mat").LEO_DATA_Apparent;

% Run Extended Kalman Filter
% Storage Matrices
XMat = []; CovMat = []; PropMat = []; ResMat = [];

% Sizes
nx = 11; nz = 2;

% Populate Filter Structures
xHatk = [IC.R_ECI_0' IC.V_ECI_0' IC.C_d 1.0 0 0 0]';      % State Estimate at k
PHatk = diag([1, 1, 1, 0.01, 0.01, 0.01, 0, 0, 0, 0, 0]); % Covariance at k
xBark = xHatk; PBark = PHatk;

% Uncertainties
Rk = zeros(3, 2, 2);
Rk(1, :, :) = diag([0.010, 0.0000005]); % Measurement Noise
Rk(2, :, :) = diag([0.005, 0.0000010]);
Rk(3, :, :) = diag([0.010, 0.0000005]);
Qk_RIC = diag([1e-12, 1e-12, 1e-12]).^2;     % Process Noise

% Iterators
curTime = 0;  % Current State Time (s)
measTime = 0; % Time of Current Measurement

% Loop Through All Measurements
for k = 1:length(MeasMat)
    fprintf("PROCESSING MEAS %d\n", k);
    curStation = MeasMat(k, 1); % Measurement Station
    measTime = MeasMat(k, 2);   % Measurement Time
    zk = MeasMat(k, 3:4)';      % Measurement (Range, Range-Rate)
    Phik = eye(nx);

    R_S_ECEFk = C.station_pos(:, curStation);
    EOP = get_EOP(IC.epoch + curTime / 86400);
    [PMat, NMat, SMat, WMat] = get_ECEF_ECI(IC.epoch + curTime / 86400, ...
                                            EOP.dUT1, EOP.xp, EOP.yp);
    PNSMat = PMat * NMat * SMat;

    % Propagate State to Measurement Time
    if (curTime < measTime)
        tVec = curTime:measTime;
        delt = measTime - curTime;
        [~, state_hist] = ode45(@state_prop_no_c_d, tVec, [xHatk; reshape(Phik, 121, 1)], ...
                                C.options, IC.epoch + curTime / 86400, getAMat, getAcc);
        Phik = reshape(state_hist(end, 12:end), 11, 11);
        xBark = state_hist(end, 1:11)';
        Q_ECI_RIC = get_ECI_RIC(xBark(1:3), xBark(4:6));
        Qk = Q_ECI_RIC' * Qk_RIC;
        PBark = Phik * PHatk * Phik' + blkdiag(delt^2 * [delt^2/4 * Qk, delt^2/2 * Qk;
                                                delt^2/2 * Qk, Qk], zeros(5));
        curTime = measTime;
    end

    % Perform Measurement Update
    HTildek = getHTildeMat(PNSMat(1, 1), PNSMat(1, 2), PNSMat(1, 3), ...
                           PNSMat(2, 1), PNSMat(2, 2), PNSMat(2, 3), ...
                           PNSMat(3, 1), PNSMat(3, 2), PNSMat(3, 3), ...
                           R_S_ECEFk(1), R_S_ECEFk(2), R_S_ECEFk(3), ... % Station Position in ECEF
                           WMat(1, 1), WMat(1, 2), WMat(1, 3), ...
                           WMat(2, 1), WMat(2, 2), WMat(2, 3), ...
                           WMat(3, 1), WMat(3, 2), WMat(3, 3), ...
                           xBark(4), xBark(5), xBark(6), ...             % Velocity
                           xBark(1), xBark(2), xBark(3));                % Position
    Gk = getGVec(PNSMat(1, 1), PNSMat(1, 2), PNSMat(1, 3), ...
                 PNSMat(2, 1), PNSMat(2, 2), PNSMat(2, 3), ...
                 PNSMat(3, 1), PNSMat(3, 2), PNSMat(3, 3), ...
                 R_S_ECEFk(1), R_S_ECEFk(2), R_S_ECEFk(3), ... % Station Position in ECEF
                 WMat(1, 1), WMat(1, 2), WMat(1, 3), ...
                 WMat(2, 1), WMat(2, 2), WMat(2, 3), ...
                 WMat(3, 1), WMat(3, 2), WMat(3, 3), ...
                 xBark(4), xBark(5), xBark(6), ...             % Velocity
                 xBark(1), xBark(2), xBark(3));                % Position
    zBark = zk - Gk; 
    KMat = PBark * HTildek' / (HTildek * PBark * HTildek' + reshape(Rk(curStation, :, :), 2, 2));
    xHatk = xBark + KMat * zBark;
    PHatk = (eye(nx) - KMat * HTildek) * PBark * (eye(nx) - KMat * HTildek)' ...
             + KMat * reshape(Rk(curStation, :, :), 2, 2) * KMat';
    
    % Store State and Covariance
    XMat = [XMat; xHatk'];
    CovMat = [CovMat; reshape(PHatk, 1, 121)];
end

fprintf("\nFinal Position: %.6f, %.6f, %.6f\n", XMat(end, 1:3));
fprintf("Final Velocity: %.6f, %.6f, %.6f\n\n", XMat(end, 4:6));

% Plot State
TimeMat = MeasMat(:, 2);
f1 = figure(1);
subplot(3, 1, 1);
scatter(TimeMat, XMat(:, 1), 'r', 'filled');
grid on
title('Estimated Position Over 24 Hours');
ylabel('x (km)');
subplot(3, 1, 2);
scatter(TimeMat, XMat(:, 2), 'r', 'filled');
grid on;
ylabel('y (km)');
subplot(3, 1, 3);
scatter(TimeMat, XMat(:, 3), 'r', 'filled');
grid on;
xlabel('Time (s)');
ylabel('z (km)');
f2 = figure(2);
subplot(3, 1, 1);
scatter(TimeMat, XMat(:, 4), 'r', 'filled');
grid on
title('Estimated Velocity Over 24 Hours');
ylabel('x (km/s)');
subplot(3, 1, 2);
scatter(TimeMat, XMat(:, 5), 'r', 'filled');
grid on;
ylabel('y (km/s)');
subplot(3, 1, 3);
scatter(TimeMat, XMat(:, 6), 'r', 'filled');
grid on;
xlabel('Time (s)');
ylabel('z (km/s)');
movegui(f1, "northwest");
movegui(f2, "southwest");

% End Timer
toc;