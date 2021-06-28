function [BCO, wx, wy, wz, T1_D, T2_D, T3_D] = EulerSequence(rotation1, rotation2, rotation3)
% EulerSequence(rotation1, rotation2, rotation3) computes the BCO matrix in
% symbolic form for a given sequence of euler angles as dictated by
% the rotation sequence of rotation 1, 2, and 3.
%
% INPUTS:
%   rotation1 = first rotation, either 'X', 'Y', 'Z', or [] for no rotation
%   rotation2 = second rotation, either 'X', 'Y', 'Z', or [] for no rotation
%   rotation3 = third rotation, either 'X', 'Y', 'Z', or [] for no rotation
% 
% OUTPUTS:
%   BCO = rotation matrix
%   wx = angular velocity in x
%   wy = angular velocity in y
%   wz = angular velocity in z
%   T1_D = theta 1 dot
%   T2_D = theta 2 dot
%   T3_D = theta 3 dot
%
% EXAMPLES:
%   Sequence YZX:
%   EulerSequence('Y', 'Z', 'X')
%   EulerSequence('YZX')
%       Results in BCO = [RX*RZ*RY] - THIS IS CORRECT!!!
%       Equal to OCB = [RY'*RZ'*RX'] - used by most people

% Modifications
% 1. CDY, 7/21/17, added handling of rotaitons will less than three
%    sequences. 
% 2. CDY, 10/18/20, standardized rotation scheme, added new "sequence"
%   input for easy. 

% handle nargin - new addition
if nargin == 1
    % parse into other rotation sequence
    if ischar(rotation1) ~= 1
        error('Must be string, single quotes.');
    end
    seq = rotation1;
    rotation1 = seq(1);
    rotation2 = seq(2);
    rotation3 = seq(3);
end

% check inputs
if isempty(rotation1) == 1
    error('rotation1 must be a non-empty matrix');
end
R1 = eye(3);
R2 = eye(3);
R3 = eye(3);
cntr = 1;

% form rotations
syms theta_1c(t) theta_2c(t) theta_3c(t)
d1c = diff(theta_1c(t), t);
d2c = diff(theta_2c(t), t);
d3c = diff(theta_3c(t), t);

% form matrices
R1 = populate(rotation1, theta_1c);
if isempty(rotation2) ~= 1
    R2 = populate(rotation2, theta_2c);
    cntr = cntr + 1;
end
if isempty(rotation3) ~= 1
    R3 = populate(rotation3, theta_3c);
    cntr = cntr + 1;
end

% calculate BCO
BCOt = R3*R2*R1;
OCB = transpose(BCOt);
OCBd = diff(OCB);
wxt = [0, 0, 1]*BCOt*OCBd*[0; 1; 0];
wyt = [1, 0, 0]*BCOt*OCBd*[0; 0; 1];
wzt = [0, 1, 0]*BCOt*OCBd*[1; 0; 0];

% remove time dependent terms
syms t_1 t_2 t_3 t_1_D t_2_D t_3_D
BCO = subs(BCOt, [theta_1c, theta_2c, theta_3c], [t_1, t_2, t_3]);
wx = simplify(subs(wxt, [theta_1c, theta_2c, theta_3c, d1c, d2c, d3c], ...
    [t_1, t_2, t_3, t_1_D, t_2_D, t_3_D]));
wy = simplify(subs(wyt, [theta_1c, theta_2c, theta_3c, d1c, d2c, d3c], ...
    [t_1, t_2, t_3, t_1_D, t_2_D, t_3_D]));
wz = simplify(subs(wzt, [theta_1c, theta_2c, theta_3c, d1c, d2c, d3c], ...
    [t_1, t_2, t_3, t_1_D, t_2_D, t_3_D]));

% extract inverse relations
syms WX WY WZ
[A, b] = equationsToMatrix([wx == WX, wy == WY, wz == WZ], [t_1_D t_2_D t_3_D]);
t1 = A\b;
T1_D = simplify(t1(1));
T2_D = simplify(t1(2));
T3_D = simplify(t1(3));

end