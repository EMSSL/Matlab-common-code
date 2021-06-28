function [t1, t2, t3, q0, q1, q2, q3] = MatrixToComponents(BCO, seq)
% returns the euler angles and quaternions
% INPUTS
%   OCB = rotation matrix
%   seq = matrix:
%           1 = 'xyz'
%           5 = 'yxz'
%          10 = 'zyx'
%          15 = 'xzy'
%
% OUTPUTS
%   t1, t2, t3      = Euler angles, rad, according to sequence
%   q0, q1, q2, q3  = quaternions 

% error handle
if isnumeric(BCO) ~= 1
    error('BCO is not a numeric matrix.');
end

% handle string inputs for sequence
seqn = seq;             % preallocate as number 
if ischar(seq) == 1
    switch seq
        case 'XYZ'
            seqn = 1;
        case 'YXZ'
            seqn = 5;
        case 'ZYX'
            seqn = 10;
        case 'XZY'
            seqn = 15;
    end
end



% do this yourself
% [BCO, wx, wy, wz, T1_D, T2_D, T3_D] = EulerSequence(rotation1, rotation2, rotation3)
switch seqn
    case 1
        % seq = 1;    % XYZ 
        % t1 = X, t2 = Y, t3 = Z
        % MakeBCO(sym('t_1'), sym('t_2'), sym('t_3'), 'XYZ')
        % EulerSequence('X', 'Y', 'Z')
        
        % BCO =
        % 
        % [  cos(t_2)*cos(t_3), cos(t_1)*sin(t_3) + cos(t_3)*sin(t_1)*sin(t_2), sin(t_1)*sin(t_3) - cos(t_1)*cos(t_3)*sin(t_2)]
        % [ -cos(t_2)*sin(t_3), cos(t_1)*cos(t_3) - sin(t_1)*sin(t_2)*sin(t_3), cos(t_3)*sin(t_1) + cos(t_1)*sin(t_2)*sin(t_3)]
        % [           sin(t_2),                             -cos(t_2)*sin(t_1),                              cos(t_1)*cos(t_2)]

        t1 = atan2(-BCO(3, 2), BCO(3, 3));      % t1 = atan2(sin(t_1), cos(t_1)) = angleX
        t2 = asin(BCO(3, 1));                   % t2 = asin(sin(t_2)) = angleY
        t3 = atan2(-BCO(2, 1), BCO(1, 1));      % t3 = atan2(sin(t_3), cos(t_3)) = angleZ
        
        % % testing the angles, should result in a zero matrix
        % t_1 = t1;
        % t_2 = t2;
        % t_3 = t3;
        % BCOT = [  cos(t_2)*cos(t_3), cos(t_1)*sin(t_3) + cos(t_3)*sin(t_1)*sin(t_2), sin(t_1)*sin(t_3) - cos(t_1)*cos(t_3)*sin(t_2);
        % -cos(t_2)*sin(t_3), cos(t_1)*cos(t_3) - sin(t_1)*sin(t_2)*sin(t_3), cos(t_3)*sin(t_1) + cos(t_1)*sin(t_2)*sin(t_3);
        %           sin(t_2),                             -cos(t_2)*sin(t_1),                              cos(t_1)*cos(t_2)];
        % BCO - BCOT
        % keyboard
        
    case 10
        % seq = 10;    % ZYX 
        % t1 = Z, t2 = Y, t3 = X
        % MakeBCO(sym('t_1'), sym('t_2'), sym('t_3'), 'ZYX')
        
        % BCO = 
        % [                             cos(t_1)*cos(t_2),                             cos(t_2)*sin(t_1),         -sin(t_2);
        % cos(t_1)*sin(t_2)*sin(t_3) - cos(t_3)*sin(t_1), cos(t_1)*cos(t_3) + sin(t_1)*sin(t_2)*sin(t_3), cos(t_2)*sin(t_3);
        % sin(t_1)*sin(t_3) + cos(t_1)*cos(t_3)*sin(t_2), cos(t_3)*sin(t_1)*sin(t_2) - cos(t_1)*sin(t_3), cos(t_2)*cos(t_3)];
        
        t1 = atan2(BCO(1, 2), BCO(1, 1));   % atan2(cos(t_2)*sin(t_1), cos(t_1)*cos(t_2)) = atan2(sin(t_1), cos(t_2))
        t2 = asin(-BCO(1, 3));
        t3 = atan2(BCO(2, 3), BCO(3, 3));   % atan2(cos(t_2)*sin(t_3), cos(t_2)*cos(t_3)) = atan2(sin(t_3), cos(t_3))
        
        % % testing the angles, should result in a zero matrix
        % t_1 = t1;
        % t_2 = t2;
        % t_3 = t3;
        % BCOT = [                             cos(t_1)*cos(t_2),                             cos(t_2)*sin(t_1),         -sin(t_2);
        %     cos(t_1)*sin(t_2)*sin(t_3) - cos(t_3)*sin(t_1), cos(t_1)*cos(t_3) + sin(t_1)*sin(t_2)*sin(t_3), cos(t_2)*sin(t_3);
        %     sin(t_1)*sin(t_3) + cos(t_1)*cos(t_3)*sin(t_2), cos(t_3)*sin(t_1)*sin(t_2) - cos(t_1)*sin(t_3), cos(t_2)*cos(t_3)];
        % BCO - BCOT
        % keyboard

        
    case 15
        % seq = 15;    % XZY
        % t1 = X, t2 = Z, t3 = Y
        % MakeBCO(sym('t_1'), sym('t_2'), sym('t_3'), 'XZY')
        
        % BCO = 
        %   [ cos(t_2)*cos(t_3), sin(t_1)*sin(t_3) + cos(t_1)*cos(t_3)*sin(t_2), cos(t_3)*sin(t_1)*sin(t_2) - cos(t_1)*sin(t_3)]
        %   [         -sin(t_2),                              cos(t_1)*cos(t_2),                              cos(t_2)*sin(t_1)]
        %   [ cos(t_2)*sin(t_3), cos(t_1)*sin(t_2)*sin(t_3) - cos(t_3)*sin(t_1), cos(t_1)*cos(t_3) + sin(t_1)*sin(t_2)*sin(t_3)]
        
        t1 = atan2(BCO(2, 3), BCO(2, 2));       % t1 = atan2(sin(t_1), cos(t_1)) = angleX
        t2 = asin(-BCO(2, 1));                  % t2 = asin(sin(t_2)) = angleZ
        t3 = atan2(BCO(3, 1), BCO(1, 1));      % t3 = atan2(sin(t_3), cos(t_3)) = angleY

        % % testing the angles, should result in a zero matrix
        % t_1 = t1;
        % t_2 = t2;
        % t_3 = t3;
        % BCOT = [ cos(t_2)*cos(t_3), sin(t_1)*sin(t_3) + cos(t_1)*cos(t_3)*sin(t_2), cos(t_3)*sin(t_1)*sin(t_2) - cos(t_1)*sin(t_3);
        %                  -sin(t_2),                              cos(t_1)*cos(t_2),                              cos(t_2)*sin(t_1);
        %          cos(t_2)*sin(t_3), cos(t_1)*sin(t_2)*sin(t_3) - cos(t_3)*sin(t_1), cos(t_1)*cos(t_3) + sin(t_1)*sin(t_2)*sin(t_3)];
        % BCO - BCOT
        % keyboard
        
    otherwise
        error('seq not recognized.');
end

% reset quats
% q0 = 0;
% q1 = 0;
% q2 = 0;
% q3 = 0;
%   BCO_q =
%       [q0c^2 + q1c^2 - q2c^2 - q3c^2,         2*q0c*q3c + 2*q1c*q2c,         2*q1c*q3c - 2*q0c*q2c]
%       [        2*q1c*q2c - 2*q0c*q3c, q0c^2 - q1c^2 + q2c^2 - q3c^2,         2*q0c*q1c + 2*q2c*q3c]
%       [        2*q0c*q2c + 2*q1c*q3c,         2*q2c*q3c - 2*q0c*q1c, q0c^2 - q1c^2 - q2c^2 + q3c^2]
% http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
% IF USING THIS LINK, HIS BCO IS MAZZOLENI'S OCB!!!
OCB = transpose(BCO);
TrOCB = OCB(1, 1) + OCB(2, 2) + OCB(3, 3);
if TrOCB > 0
    % normal, all good
    % S = 0.5/sqrt(TrOCB + 1);
    % S = 0.5/sqrt(TrOCB);
    S = 2*sqrt(TrOCB + 1);
    q0 = 0.25*S;
    q1 = (OCB(3, 2) - OCB(2, 3))/S;
    q2 = (OCB(1, 3) - OCB(3, 1))/S;
    q3 = (OCB(2, 1) - OCB(1, 2))/S;
elseif (OCB(1, 1) > OCB(2, 2)) && (OCB(1, 1) > OCB(3, 3))
    % (1, 1) is biggest
    S = 2*sqrt(1 + OCB(1, 1) - OCB(2, 2) - OCB(3, 3));
    q0 = (OCB(3, 2) - OCB(2, 3))/S;
    q1 = 0.25*S;
    q2 = (OCB(1, 2) + OCB(2, 1))/S;
    q3 = (OCB(1, 3) + OCB(3, 1))/S;
elseif OCB(2, 2) > OCB(3, 3)
    % (2, 2) is biggest
    S = 2*sqrt(1 + OCB(2, 2) - OCB(1, 1) - OCB(3, 3));
    q0 = (OCB(1, 3) - OCB(3, 1))/S;
    q1 = (OCB(1, 2) + OCB(2, 1))/S;
    q2 = 0.25*S;
    q3 = (OCB(2, 3) + OCB(3, 2))/S;
else 
    % (3, 3) is biggest
    S = 2*sqrt(1 + OCB(3, 3) - OCB(1, 1) - OCB(2, 2));
    q0 = (OCB(2, 1) - OCB(1, 2))/S;
    q1 = (OCB(1, 3) + OCB(3, 1))/S;
    q2 = (OCB(2, 3) + OCB(3, 2))/S;
    q3 = 0.25*S;
end

% % test with matlab functions
% [t1, t2, t3]
% eul = rotm2eul(BCO', "ZYX")
% keyboard

end       
        