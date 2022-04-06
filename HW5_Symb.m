%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASE 389.4 Methods of Orbit Determination
% Homework 5
% Author: Alex Chiu (ac68767)
% Last Edited: 04/04/2022
% Summary: Symbolic Code for HW5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define A Matrix
syms A [11 11]

% Define State Variables
syms x y z vx vy vz C_D C_solar b1 b2 b3

% Radius Function
syms r(x, y, z)
r(x, y, z) = sqrt(x^2 + y^2 + z^2);

% Define Acceleration
syms A_Total [3 1]
syms ax ay az

% Combine Values into Vector 
syms VecForm(v1, v2, v3);
VecForm(v1, v2, v3) = [v1; v2; v3];

% Two Body Acceleration
syms mu_e
syms ax_2b ay_2b az_2b
syms Acc_2B [3 1]
ax_2b = -mu_e / r(x, y, z)^3 * x;
ay_2b = -mu_e / r(x, y, z)^3 * y;
az_2b = -mu_e / r(x, y, z)^3 * z;
Acc_2B = VecForm(ax_2b, ay_2b, az_2b);
A_Total = Acc_2B;

% J2 Acceleration
syms J2 R_e
syms ax_J2 ay_J2 az_J2
ax_J2 = 3*mu_e*J2 * R_e^2/(2*r(x, y, z)^5) * x*(5*(z/r(x, y, z))^2 - 1);
ay_J2 = 3*mu_e*J2 * R_e^2/(2*r(x, y, z)^5) * y*(5*(z/r(x, y, z))^2 - 1);
az_J2 = 3*mu_e*J2 * R_e^2/(2*r(x, y, z)^5) * z*(5*(z/r(x, y, z))^2 - 3);
Acc_J2 = VecForm(ax_J2, ay_J2, az_J2);
A_Total = A_Total + Acc_J2;

% Drag Acceleration
syms A_sp m rho_0 r_0 H w_e
syms ax_d ay_d az_d
syms v_rel_norm(x, y, vx, vy, vz)
v_rel_norm(x, y, vx, vy, vz) = sqrt((vx+w_e*y)^2 + (vy-w_e*x)^2 + vz^2); % TODO
ax_d = -1/2*C_D*A_sp/m*rho_0*exp(-(r(x, y, z)-r_0)/H)*1000*v_rel_norm(x, y, vx, vy, vz)*(vx+w_e*y);
ay_d = -1/2*C_D*A_sp/m*rho_0*exp(-(r(x, y, z)-r_0)/H)*1000*v_rel_norm(x, y, vx, vy, vz)*(vy-w_e*x);
az_d = -1/2*C_D*A_sp/m*rho_0*exp(-(r(x, y, z)-r_0)/H)*1000*v_rel_norm(x, y, vx, vy, vz)*vz;
Acc_D = VecForm(ax_d, ay_d, az_d);
A_Total = A_Total + Acc_D;

% Populate Vectors
X = [x y z vx vy vz C_D C_solar b1 b2 b3];
F = [vx vy vz A_Total(1) A_Total(2) A_Total(3) 0 0 0 0 0];

% Find Partial Derivatives
for i = 1:11
    for j = 1:11
        A(i, j) = diff(F(i), X(j));
    end
end

% Define H_Tilde
syms H_Tilde [2 11]
syms R_S_ECEF R_S_ECI [3 1]
syms range(x, y, z) range_rate(x, y, z, vx, vy, vz)
W = sym('W', [3 3]);
PNS = sym('PNS', [3 3]);

R_S_PEF = W * R_S_ECEF;
R_S_ECI = PNS * R_S_PEF;
V_S_ECI = PNS * (cross([0; 0; w_e], R_S_PEF));

range(x, y, z) = sqrt((x - R_S_ECI(1))^2 + (y - R_S_ECI(2))^2 + (z - R_S_ECI(3))^2);
range_rate(x, y, z, vx, vy, vz) = ((x - R_S_ECI(1))*(vx - V_S_ECI(1)) + ...
                                   (y - R_S_ECI(2))*(vy - V_S_ECI(2)) + ...
                                   (z - R_S_ECI(3))*(vz - V_S_ECI(3))) / range(x, y, z);

Z = [range(x, y, z); range_rate(x, y, z, vx, vy, vz)];

% Find Partial Derivatives
for i = 1:2
    for j = 1:11
        H_Tilde(i, j) = diff(Z(i), X(j));
    end
end

% Load Constants
mu_e = C.mu_e;
J2 = C.J2;
R_e = C.R_e;
A_sp = C.A_sp;
m = C.m;
rho_0 = C.rho0;
r_0 = C.r0;
H = C.H;
w_e = C.w_e;

% Sub Values
Acc_2B = subs(Acc_2B);
Acc_J2 = subs(Acc_J2);
Acc_D = subs(Acc_D);
A = subs(A);
F = subs(F);
H_Tilde = subs(H_Tilde);
Z = subs(Z);

% Define Final Functions
getAcc2B = matlabFunction(Acc_2B);
getAccJ2 = matlabFunction(Acc_J2);
getAccD = matlabFunction(Acc_D);
getAMat = matlabFunction(A);
getF = matlabFunction(F);
getHMat = matlabFunction(H_Tilde);
getZ = matlabFunction(Z);