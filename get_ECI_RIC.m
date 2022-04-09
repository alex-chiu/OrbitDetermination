function [Q_ECI_RIC] = get_ECI_RIC(R_ECI, V_ECI)
% get_ECI_Body : Get transformation matrix from ECI to RIC frame
%
% INPUTS
%
% R_ECI ------ 3-by-1 position of spacecraft in ECI
%
% V_ECI ------ 3-by-1 velocity of spacecraft in ECI
%
% OUTPUTS
%
% Q_ECI_RIC - Rotation matrix from Body to ECI frame.
% 
%+------------------------------------------------------------------------------+
% References: Lecture Notes 
%
% Author: Alex Chiu
%
% Last Edited: 4/9/2022
%+==============================================================================+

Radial = R_ECI / norm(R_ECI);
Intrack = V_ECI / norm(V_ECI);
Crosstrack = cross(Radial, Intrack);

Q_ECI_RIC = [Radial Intrack Crosstrack]';
