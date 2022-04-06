function [R] = eulerRotate(axis, angle)
% eulerRotate : Converts a single Euler angle (in radians) into a direction cosine matrix.
%
% INPUTS
%
% axis ------- axis about which to perform the rotation (1, 2, or 3)
%
% angle ------ angle of rotation
%
% OUTPUTS
%
% R ---------- 3-by-3 direction cosine matrix 
% 
%+------------------------------------------------------------------------------+
% References: Lecture Notes 
%
% Author: Alex Chiu
%
% Last Edited: 1/21/2022
%+==============================================================================+

if (axis == 1)
    R = [1, 0, 0;
         0, cos(angle), sin(angle);
         0, -sin(angle), cos(angle)];
elseif (axis == 2)
    R = [cos(angle), 0, -sin(angle);
         0, 1, 0;
         sin(angle), 0, cos(angle)];
elseif (axis == 3)
    R = [cos(angle), sin(angle), 0;
         -sin(angle), cos(angle), 0;
         0, 0, 1];
else
    error('Axis has to be 1, 2, or 3');
end

end