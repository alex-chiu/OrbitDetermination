function [EOP] = get_EOP(JD)
% get_EOP : Get Earth Orientation Parameters for ECEF-ECI transformation
%
% INPUTS
%
% JD-- ------ Time of interest, as a Julian Date in UTC
%
% OUTPUTS
%
% EOP ------- Structure containing all EOPs
% 
%+------------------------------------------------------------------------------+
% References: Lecture Notes 
%
% Author: Alex Chiu
%
% Last Edited: 4/8/2022
%+==============================================================================+

B361 = readmatrix(".\Data\bulletinb-361.csv");
B362 = readmatrix(".\Data\bulletinb-361.csv");
B363 = readmatrix(".\Data\bulletinb-363.csv");
B364 = readmatrix(".\Data\bulletinb-364.csv");
data = [B361; B362; B363; B364];

MJD = JD - 2400000.5;

i = 1;
while (MJD > data(i, 1))
    i = i + 1;
end

EOP.dUT1 = (data(i - 1, 11) + ((MJD - data(i - 1, 1))/(data(i, 1) - data(i - 1, 1)) ...
           * (data(i, 11) - data(i - 1, 11)))) / 1000;
EOP.xp = (data(i - 1, 6) + ((MJD - data(i - 1, 1))/(data(i, 1) - data(i - 1, 1)) ...
           * (data(i, 6) - data(i - 1, 6)))) / 1000;
EOP.yp = (data(i - 1, 8) + ((MJD - data(i - 1, 1))/(data(i, 1) - data(i - 1, 1)) ...
           * (data(i, 8) - data(i - 1, 8)))) / 1000;