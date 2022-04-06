%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASE 389.4 Methods of Orbit Determination
% Homework 5
% Author: Alex Chiu (ac68767)
% Last Edited: 04/04/2022
% Summary: Main Computation Code for HW5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
close all;
format long g;

tic

% Load Constants
constantsScript;

% Problem 1
% Derive A and H_Tilde Matrices 

R0 = [6990077.798814194 1617465.311978378 22679.810569245355]' / 1000; % km
V0 = [-1675.13972506056 7273.72441330686 252.688512916741]' / 1000; % km / s

% Load Symbolic Matrix
HW5_Symb;

% Plug Values In
A0 = getAMat(C.C_d, V0(1), V0(2), V0(3), R0(1), R0(2), R0(3));
R_S_ECEF = C.station_pos(:, 1);
[PMat, NMat, SMat, WMat] = getECEF2ECI(2458150.70833, 0.1965014, 0.015361, 0.288259);
W = WMat;
PNS = PMat * NMat * SMat;

H0 = getHMat(PNS(1, 1), PNS(1, 2), PNS(1, 3), ...
             PNS(2, 1), PNS(2, 2), PNS(2, 3), ...
             PNS(3, 1), PNS(3, 2), PNS(3, 3), ...
             R_S_ECEF(1), R_S_ECEF(2), R_S_ECEF(3), ...
             W(1, 1), W(1, 2), W(1, 3), ...
             W(2, 1), W(2, 2), W(2, 3), ...
             W(3, 1), W(3, 2), W(3, 3), ...
             V0(1), V0(2), V0(3), R0(1), R0(2), R0(3));

% Compare A Matrices
A_sol = load("Data\A_t0.mat").A;
A_diff = abs((A0(1:7, 1:7) - A_sol) ./ A_sol);

% Compare H_Tilde Matrices
H_sol = load("Data\H_Tilde_t0.mat").H_TILDA;
H_diff = abs((H0(:, 1:7) - H_sol) ./ H_sol);

disp(H_diff);

% Problem 2
Phi0 = eye(7);
X0 = [R0' V0' reshape(Phi0(:), 1, 49)]';
t = 0:60:21600;

% Set ODE Options
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

[t_hist, y_hist] = ode45(@full_prop, t, X0, options, 0, getAcc2B, getAccJ2, getAccD, getAMat, C);

% Compare Phi Matrices
Phi_sol = load("Data\Phi_21600_0.mat").PHI_t_120;
Phi_diff = abs((reshape(y_hist(end, 7:end), 7, 7) - Phi_sol) ./ Phi_sol);

truth = load("Data\lighttime_truth.mat").lighttime_truth;

% Plot All Figures
f1 = figure(1);
histogram(reshape(log10(abs(double(A_diff))), 49, 1), 4);
grid on;
title("Rel Diff of A at Initial Time");

f2 = figure(2);
subplot(3, 1, 1);
plot(t_hist, y_hist(:, 1));
grid on;
hold on;
plot(truth(:, 7), truth(:, 1));
ylabel("x (km)");
title("Position Error History");
legend("Propagated", "Truth");
subplot(3, 1, 2);
plot(t_hist, y_hist(:, 2));
grid on;
hold on;
plot(truth(:, 7), truth(:, 2));
ylabel("y (km)");
subplot(3, 1, 3);
plot(t_hist, y_hist(:, 3));
grid on;
hold on;
plot(truth(:, 7), truth(:, 3));
xlabel("Time (s)");
ylabel("z (km)");

f3 = figure(3);
histogram(reshape(log10(abs(double(Phi_diff))), 49, 1), 4);
grid on;
title("Rel Diff of STM at Final Time");

movegui(f1, "northwest");
movegui(f2, "southwest");
movegui(f3, "north");

toc