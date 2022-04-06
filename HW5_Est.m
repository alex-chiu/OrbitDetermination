%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASE 389.4 Methods of Orbit Determination
% Homework 5
% Author: Alex Chiu (ac68767)
% Last Edited: 04/05/2022
% Summary: Main Estimation Code for HW5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
close all;
format long g;

tic

% Load Constants and Data
constantsScript;
MeasMat = load("Data\LEO_DATA_Apparent.mat").LEO_DATA_Apparent;

dtIn = 60;
tVec = 0:dtIn:21600; % 6 Hours

% Setup
N = length(tVec);
HW5_Symb;
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

% Initial Values
R0 = [6990077.798814194 1617465.311978378 22679.810569245355]' / 1000; % km
V0 = [-1675.13972506056 7273.72441330686 252.688512916741]' / 1000; % km / s

XMat = [];
PropMat = [];
PropTimeMat = [];
CovMat = [];
TimeMat = [];
ResMat = [];

xHatk = [R0' V0' 1.88 1.0 0 0 0]';
xBarkp1 = xHatk;
Phik = eye(11);
Phikp1 = Phik;
Pk = diag([500, 500, 500, 50, 50, 50, 1, 1, 1, 1, 1]);
PBarkp1 = Pk;
Rk = diag([0.005, 0.000001]);

Qk = diag(1e-21 * ones(6, 1));

RangeErr = 0;
RangeRateErr = 0;

curTime = 0;
measIdx = 1;
while curTime < tVec(end)
    nextMeasTime = MeasMat(measIdx, 2);
    
    % Propagate to Next Measurement
    if (curTime < nextMeasTime)
        % fprintf("PROPAGATING FROM %d TO %d\n", curTime, nextMeasTime);
        propTimeVec = curTime:nextMeasTime;
        delt = nextMeasTime - curTime;
        [t_hist, y_hist] = ode45(@state_prop, propTimeVec, [xHatk; reshape(Phik, 121, 1)], options, 0, getAcc2B, getAccJ2, getAccD, getAMat);
        PropMat = [PropMat; y_hist(1:end, 1:11)];
        PropTimeMat = [PropTimeMat; t_hist(1:end)];
        xBarkp1 = y_hist(end, 1:11)';
        Phikp1 = reshape(y_hist(end, 12:end), 11, 11);
        Ak = getAMat(xHatk(7), xHatk(4), xHatk(5), xHatk(6), xHatk(1), xHatk(2), xHatk(3));
        Gammak = zeros(6, 6);
        Gammak(1:3, 1:3) = diag([delt/2, delt/2, delt/2]);
        Gammak(4:6, 1:3) = diag(ones(3, 1));
        Gammak(1:3, 4:6) = diag(ones(3, 1));
        PBarkp1 = Ak * Pk * Ak' + delt^2 * blkdiag(Gammak*Qk*Gammak', zeros(5));
        curTime = nextMeasTime;
    else
        zkp1 = MeasMat(measIdx, 3:4)';
        R_S_ECEF = C.station_pos(:, MeasMat(measIdx, 1));
        [PMat, NMat, SMat, WMat] = getECEF2ECI(C.start_time+curTime/86400, 0.1965014, 0.015361, 0.288259);
        W = WMat;
        PNS = PMat * NMat * SMat;
        zBarkp1 = getZ(PNS(1, 1), PNS(1, 2), PNS(1, 3), ...
                 PNS(2, 1), PNS(2, 2), PNS(2, 3), ...
                 PNS(3, 1), PNS(3, 2), PNS(3, 3), ...
                 R_S_ECEF(1), R_S_ECEF(2), R_S_ECEF(3), ...
                 W(1, 1), W(1, 2), W(1, 3), ...
                 W(2, 1), W(2, 2), W(2, 3), ...
                 W(3, 1), W(3, 2), W(3, 3), ...
                 xBarkp1(4), xBarkp1(5), xBarkp1(6), xBarkp1(1), xBarkp1(2), xBarkp1(3));
        nukp1 = zkp1 - zBarkp1;
        ResMat = [ResMat; zBarkp1'];
        RangeErr = RangeErr + nukp1(1)^2;
        RangeRateErr = RangeRateErr + nukp1(2)^2;
        % fprintf("MEASUREMENT DIFFERENCE: %f, %f\n", nukp1);
        
        Hkp1 = getHMat(PNS(1, 1), PNS(1, 2), PNS(1, 3), ...
                       PNS(2, 1), PNS(2, 2), PNS(2, 3), ...
                       PNS(3, 1), PNS(3, 2), PNS(3, 3), ...
                       R_S_ECEF(1), R_S_ECEF(2), R_S_ECEF(3), ...
                       W(1, 1), W(1, 2), W(1, 3), ...
                       W(2, 1), W(2, 2), W(2, 3), ...
                       W(3, 1), W(3, 2), W(3, 3), ...
                       xBarkp1(4), xBarkp1(5), xBarkp1(6), xBarkp1(1), xBarkp1(2), xBarkp1(3));
        Skp1 = Hkp1 * PBarkp1 * Hkp1' + Rk;
        Wkp1 = PBarkp1 * Hkp1' / Skp1;
        xHatkp1 = xBarkp1 + Wkp1 * nukp1;
        Pkp1 = PBarkp1 - Wkp1 * Skp1 * Wkp1';
        
        XMat = [XMat; xHatk'];
        CovMat = [CovMat; reshape(Pk, 1, 121)];
        TimeMat = [TimeMat; curTime];
    
        xHatk = xHatkp1;
        Phik = Phikp1;
        Pk = Pkp1;
    
        measIdx = measIdx + 1;
    end
end

RMS_R = sqrt(RangeErr / length(MeasMat));
RMS_RR = sqrt(RangeRateErr / length(MeasMat));

fprintf("\nRMS for Range: %f\n", RMS_R);
fprintf("RMS for Range Rate: %f\n", RMS_RR);

% Post Fit Residuals

f1 = figure(1);
subplot(3, 1, 1);
scatter(TimeMat, XMat(:, 1), 'r', 'filled');
grid on
hold on
scatter(PropTimeMat, PropMat(:, 1), 1, 'b');
title('Estimated Position Over 6 Hours');
ylabel('x (km)');
subplot(3, 1, 2);
scatter(TimeMat, XMat(:, 2), 'r', 'filled');
grid on;
hold on
scatter(PropTimeMat, PropMat(:, 2), 1, 'b');
ylabel('y (km)');
subplot(3, 1, 3);
scatter(TimeMat, XMat(:, 3), 'r', 'filled');
grid on;
hold on
scatter(PropTimeMat, PropMat(:, 3), 1, 'b');
xlabel('Time (s)');
ylabel('z (km)');

f2 = figure(2);
subplot(3, 1, 1);
scatter(TimeMat, XMat(:, 4), 'r', 'filled');
grid on
hold on
scatter(PropTimeMat, PropMat(:, 4), 1, 'b');
title('Estimated Velocity Over 6 Hours');
ylabel('x (km/s)');
subplot(3, 1, 2);
scatter(TimeMat, XMat(:, 5), 'r', 'filled');
grid on;
hold on
scatter(PropTimeMat, PropMat(:, 5), 1, 'b');
ylabel('y (km/s)');
subplot(3, 1, 3);
scatter(TimeMat, XMat(:, 6), 'r', 'filled');
grid on;
hold on
scatter(PropTimeMat, PropMat(:, 6), 1, 'b');
xlabel('Time (s)');
ylabel('z (km/s)');

f3 = figure(3);
subplot(2, 1, 1);
scatter(TimeMat, MeasMat(1:end-1, 3) - ResMat(:, 1), 5, 'filled');
grid on;
ylabel('Range (km)');
title("Post Fit Residuals");
subplot(2, 1, 2);
scatter(TimeMat, MeasMat(1:end-1, 4) - ResMat(:, 2), 5, 'filled');
grid on;
xlabel('Time (s)');
ylabel('Range Rate (km/s)');

movegui(f1, "northwest");
movegui(f2, "southwest");
movegui(f3, "north");