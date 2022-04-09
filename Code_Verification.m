%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASE 389.4 Methods of Orbit Determination
% Final Project
% Author: Alex Chiu (ac68767)
% Last Edited: 04/09/2022
% Summary: Script for Testing Code Performance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear Workspace
clc; clear all; close all; format long g;

% Start Timer
tic;

% Load Constants and Initial Conditions
constants_script; initial_conditions;

% Get Function Handles for Various Matrices
[getAMat, getHTildeMat, getGVec, getAcc] = get_symb_func(C, 11);

% Flag for testing initial conditions
test_init = false;
if (test_init)
    EOP = get_EOP(IC_Test.epoch);
    [PMat, NMat, SMat, WMat] = get_ECEF_ECI(IC_Test.epoch, EOP.dUT1, EOP.xp, EOP.yp); %#ok<*UNRCH> 
    PNSMat = PMat * NMat * SMat;
    R_S_ECEF = C.station_pos(:, 1);
    
    A_0 = getAMat(IC_Test.C_d, ...                                                % C_d
                  IC_Test.V_ECI_0(1), IC_Test.V_ECI_0(2), IC_Test.V_ECI_0(3), ... % Velocity
                  IC_Test.R_ECI_0(1), IC_Test.R_ECI_0(2), IC_Test.R_ECI_0(3));    % Position
    
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
                       IC_Test.V_ECI_0(1), IC_Test.V_ECI_0(2), IC_Test.V_ECI_0(3), ... % Velocity
                       IC_Test.R_ECI_0(1), IC_Test.R_ECI_0(2), IC_Test.R_ECI_0(3));    % Position
    
    H_sol = load("..\HW 5\Data\H_Tilde_t0.mat").H_TILDA;
    H_diff = abs((H_0(:, 1:7) - H_sol) ./ H_sol);
    disp(H_diff);
end

% Load Measurements
MeasMat = load("..\HW 5\Data\LEO_DATA_Apparent.mat").LEO_DATA_Apparent;

% Run Extended Kalman Filter
% Storage MatrIC_Testes
XMat = []; CovMat = []; PropMat = []; ResMat = [];

% Sizes
nx = 11; nz = 2;

% Populate Filter Structures
xHatk = [IC_Test.R_ECI_0' IC_Test.V_ECI_0' IC_Test.C_d 1.0 0 0 0]';      % State Estimate at k
PHatk = diag([1, 1, 1, 0.01, 0.01, 0.01, 0.1, 0, 0, 0, 0]); % Covariance at k
xBark = xHatk; PBark = PHatk;

% Uncertainties
Rk = diag([0.005, 0.000001]);                       % Measurement Noise
Qk = diag([1e-12, 1e-9, 1e-12]).^2; % Process Noise

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

    % Propagate State to Measurement Time
    if (curTime < measTime)
        tVec = curTime:measTime;
        delt = measTime - curTime;
        [~, state_hist] = ode45(@state_prop, tVec, [xHatk; reshape(Phik, 121, 1)], C.options, 0, getAMat, getAcc);
        Phik = reshape(state_hist(end, 12:end), 11, 11);
        xBark = state_hist(end, 1:11)';
        PBark = Phik * PHatk * Phik' + blkdiag(delt^2 * [delt^2/4 * Qk, delt^2/2 * Qk;
                                               delt^2/2 * Qk, Qk], zeros(5));
        curTime = measTime;
    end

    % Perform Measurement Update
    R_S_ECEFk = C.station_pos(:, curStation);
    EOP = get_EOP(IC_Test.epoch + curTime / 86400);
    [PMat, NMat, SMat, WMat] = get_ECEF_ECI(IC_Test.epoch + curTime / 86400, ...
                                            EOP.dUT1, EOP.xp, EOP.yp);
    PNSMat = PMat * NMat * SMat;
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
    KMat = PBark * HTildek' / (HTildek * PBark * HTildek' + Rk);
    xHatk = xBark + KMat * zBark;
    PHatk = (eye(nx) - KMat * HTildek) * PBark * (eye(nx) - KMat * HTildek)' ...
             + KMat * Rk * KMat';
    
    % Store State and Covariance
    XMat = [XMat; xHatk'];
    CovMat = [CovMat; reshape(PHatk, 1, 121)];
end

fprintf("Final Position: %.6f, %.6f, %.6f\n", XMat(end, 1:3));
fprintf("Final Velocity: %.6f, %.6f, %.6f\n", XMat(end, 4:6));

% Plot State
TimeMat = MeasMat(:, 2);
f1 = figure(1);
subplot(3, 1, 1);
scatter(TimeMat, XMat(:, 1), 'r', 'filled');
grid on
title('Estimated Position Over 6 Hours');
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
title('Estimated Velocity Over 6 Hours');
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