function rotMat = populate(axis, angle)
% outs = populate(axis, angle) populates a rotation matrix about axis with
% angle.
%
% INPUT
%   axis = axis, either 'X', 'Y', 'Z'
%   angle = symbolic (or numeric given in radians) angle value 
% 
% OUTPUT
%   rotMat = rotation matrix of BCO (Mazzoleni notation)
%
% Written by Christopher D. Yoder on 5/30/17

% inputs
if ischar(axis) ~= 1
    error('Axis must be a string.');
end

% case for angles and matrices
if isa(angle, 'symfun') == 1 || isa(angle, 'sym') == 1
    rotMat = sym('rotMat', [3, 3]);
end
switch axis
    
    case 'X'
        rotMat(1, 1) = 1;
        rotMat(1, 2) = 0;
        rotMat(1, 3) = 0;
        rotMat(2, 1) = 0;
        rotMat(2, 2) = cos(angle);
        rotMat(2, 3) = sin(angle);
        rotMat(3, 1) = 0;
        rotMat(3, 2) = -sin(angle);
        rotMat(3, 3) = cos(angle);
        
    case 'Y'
        rotMat(1, 1) = cos(angle);
        rotMat(1, 2) = 0;
        rotMat(1, 3) = -sin(angle);
        rotMat(2, 1) = 0;
        rotMat(2, 2) = 1;
        rotMat(2, 3) = 0;
        rotMat(3, 1) = sin(angle);
        rotMat(3, 2) = 0;
        rotMat(3, 3) = cos(angle);
        
    case 'Z'
        rotMat(1, 1) = cos(angle);
        rotMat(1, 2) = sin(angle);
        rotMat(1, 3) = 0;
        rotMat(2, 1) = -sin(angle);
        rotMat(2, 2) = cos(angle);
        rotMat(2, 3) = 0;
        rotMat(3, 1) = 0;
        rotMat(3, 2) = 0;
        rotMat(3, 3) = 1;
        
    otherwise
        error('Undefined input to axis');
end

end