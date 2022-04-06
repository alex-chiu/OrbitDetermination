function [PMat, NMat, SMat, WMat] = getECEF2ECI(JD_UTC, dUT1, xp, yp)
% getECEF2ECI : Get ECEF2ECI (ITRF -> J2000) Matrix
%
% INPUTS
%
% JD_UTC ----- Julian Date in UTC
%
% dUT1 ------- Difference of UT1 from UTC in seconds
%
% xp --------- EOP for x in arcseconds
%
% yp --------- EOP for y in arcseconds
%
% OUTPUTS
%
% PMat, NMat, SMat, WMat - Various conversion matrices
%
%+------------------------------------------------------------------------------+
% References: Lecture Notes 
% https://hpiers.obspm.fr/eop-pc/index.php
%
% Author: Alex Chiu
%
% Last Edited: 4/3/2022
%+==============================================================================+   

% Constants
sec2rad = 4.848136811095359935899141e-6;

% Input Conditions
dAT = 37; % Sec
JD_TAI = JD_UTC + dAT / 86400;
JD_TT = JD_TAI + 32.184 / 86400;

% Polar Motion
xp = xp * sec2rad;
yp = yp * sec2rad;
WMat = getPolarMotion(xp, yp);

% Nutation
T_TT = (JD_TT - 2451545.0) / 36525;
ddPsi = 0; ddEps = 0;
[NMat, FA, dPsi, ~, meanEps, ~] = getNutationIAU76(T_TT, 'nut80.dat', ddPsi, ddEps);

% Sidereal Time
[SMat, ~, ~] = getSiderealTime(JD_UTC, dUT1, dPsi, meanEps, FA);

% Precession
[PMat] = getPrecessionIAU76(T_TT);

end

function [PMat] = getPrecessionIAU76(T_TT)
% INPUTS
%
% T_TT ------- Julian centuries in Terrestrial Time (TT) elapsed from
%              J2000.0
%
% OUTPUTS
%
% PMat ------- 3-by-3 precession rotation matrix 

% Constants
sec2rad = 4.848136811095359935899141e-6;
PA = [2306.2181  0.30188  0.017998;
      2004.3109 -0.42665 -0.041833;
      2306.2181  1.09468  0.018203]; % Arcseconds

PA = PA * [T_TT T_TT^2 T_TT^3]' * sec2rad; % Radians
PMat = eulerRotate(3, PA(1)) * eulerRotate(2, -PA(2)) * eulerRotate(3, PA(3));
PMat = orthodcm(PMat);
end

function [NMat, FA, dPsi, dEps, meanEps, trueEps] = getNutationIAU76(T_TT, filepath, ddPsi, ddEps)
% INPUTS
%
% T_TT ------- Julian centuries in Terrestrial Time (TT) elapsed from
%              J2000.0
%
% filepath --- String containing filepath to Nutation data file, 106 x 9 
%              set of parameters
%
% ddPsi ------ Correction to dPsi term (rad)
%
% ddEps ------ Correction to dEps term (rad)
%
% OUTPUTS
%
% NMat ------- 3-by-3 nutation rotation matrix
%
% FA --------- 5-by-1 vector of fundamental arguments
%
% dPsi ------- Nutation in longitude (rad)
%
% dEps ------- Nutation in obliquity (rad)
%
% meanEps ---- Mean obliquity of the ecliptic (rad)
%
% trueEps ---- True obliquity of the ecliptic (rad)

% Constants
sec2rad = 4.848136811095359935899141e-6;

DA = [485868.249036 1717915923.2178 31.8792 0.051635 -0.00024470;
      1287104.79305 129596581.0481 -0.5532  0.000136 -0.00001149;
      335779.526232 1739527262.8478 -12.7512 -0.001037 0.00000417;
      1072260.70369 1602961601.2090 -6.3706 0.006593 -0.00003169;
      450160.398036 -6962890.5431 7.4722 0.007702 -0.00005939];

% Load Data
NData = table2array(readtable(filepath));

% Find Mean Obliquity
meanEps = (84381.448 - 46.8150*T_TT - 0.00059*T_TT^2 + 0.001813*T_TT^3) * sec2rad;

% Find Fundamental Arguments
FA = DA * [1 T_TT T_TT^2 T_TT^3 T_TT^4]' * sec2rad; % Radians
FA = mod(FA, 2*pi);

% Find dPsi and dEps
dPsi = ddPsi;
dEps = ddEps;
for i = 1:106
    phi = NData(i, 1:5) * FA; 
    dPsi = dPsi + (NData(i, 6) + NData(i, 7) * T_TT) * sin(phi) / 10000 * sec2rad;
    dEps = dEps + (NData(i, 8) + NData(i, 9) * T_TT) * cos(phi) / 10000 * sec2rad;
end

% Find True Obliquity
trueEps = meanEps + dEps;

% Find Nutation Rotation Matrix
NMat = eulerRotate(1, -meanEps) * eulerRotate(3, dPsi) * eulerRotate(1, trueEps);
NMat = orthodcm(NMat);
end

function [WMat] = getPolarMotion(xp, yp)
% INPUTS
%
% xp --------- EOP in x direction (rad)
%
% yp --------- EOP in y direction (rad)
%
% OUTPUTS
%
% WMat ------- 3-by-3 polar motion rotation matrix 

% Find Polar Motion Rotation Matrix
WMat = eulerRotate(1, yp) * eulerRotate(2, xp);
WMat = orthodcm(WMat);
end

function [SMat, GMST, GAST] = getSiderealTime(JD_UTC, dUT1, dPsi, meanEps, FA)
% INPUTS
%
% JD_UTC ----- Julian date in UTC (days)
%
% dUT1 ------- Number of seconds between UTC and UT1 (sec)
%
% dPsi ------- Nutation in longitude (rad)
%
% FA --------- 5-by-1 vector of fundamental arguments
%
% meanEps ---- Mean obliquity of the ecliptic (rad)
%
% OUTPUTS
%
% SMat ------- 3-by-3 sidereal time rotation matrix 

% Constants
sec2rad = 4.848136811095359935899141e-6;

% Find UT1
JD_UT1 = JD_UTC + dUT1 / 86400;
T_UT1 = (JD_UT1 - 2451545.0) / 36525; % Julian Centuries
% Find Greenwich Mean Sidereal Time
GMST = 67310.54841 + (876600 * 3600 + 8640184.812866) * T_UT1 + ...
       0.093104 * T_UT1^2 - 6.2e-6 * T_UT1^3; % Seconds
GMST = mod(GMST, 86400) / 240 / 180 * pi; % Radians

% Equation of Equinoxes
Eq = dPsi * cos(meanEps) + (0.00264 * sin(FA(5)) + 0.000063 * sin(2*FA(5))) * sec2rad;

% Greenwich Apparent Sidereal Time
GAST = GMST + Eq;

SMat = eulerRotate(3, -GAST);
SMat = orthodcm(SMat);
end