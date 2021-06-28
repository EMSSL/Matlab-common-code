% Eraser_derivation_rev3
% Derivation of a flying eraser to show 6DOF example
% CDY 
%
% rev table
% rev0 - 01/17/19 - initial
% rev1 - 06/14/21 - revised for moon eraser problem
% rev2 - 06/26/21 - added missing populate fuction for S. Agrawal
% rev3 - 06/27/21 - modified to quaternions instead of euler angles 

% TO DO
% 1. DONE! - change strcmp to isa for symbolic v numeric check
% 2.       - fix handle of symbolic for eulers 
% 3.       - trace the origin of the quat to euler function expressions
% 4.       - write example code to rotate using both euler and quats, compare,
%           confirm conversion between the two. 
% 5.       - Confirm MatrixToEuler function 



clc
close all
clear
eomsfile = 'eraser_equations_rev3.txt';

% Problem Formulation
% 
% inertial frame O
%   iO = downrange, forward
%   jO = crossrange, left 
%   kO = altitude, Up
% relative frame B
%   iB = along Lx
%   jB = along Ly
%   kB = along Lz
%
% Notation:
%             r_BO_Q = linear position of Point B wrt Point O expressed in Q frame 
%           O_v_BO_B = linear velocity of B wrt O expressed in B frame
%           O_a_BO_B = linear acceleration of B wrt O expressed in B frame
%             w_OB_B = angular velocity of O frame wrt B frame expressed in B frame
%             a_OB_B = angular acceleration of O frame wrt B frame expressed in B frame
%              I_B_B = inertia tensor about Point B expressed in B frame
%      O_h_O_B_sys_B = inertial angular momentum of the system with a velocity referenced to Point O about Point B expressed in the B frame 
% Oddt_O_h_O_B_sys_B = inertial derivative of the inertial angular momentum of the system with a velocity referenced to Point O about Point B expressed in the B frame 

% derive parameters
syms t d g m Lx Ly Lz
syms xBc(t) yBc(t) zBc(t) xB yB zB              % states of point B on eraser (taken as B = CM)
syms exc(t) eyc(t) ezc(t) ex ey ez              % euler angles 
syms wxc(t) wyc(t) wzc(t) wx wy wz              % angular velocities 
syms uBc(t) vBc(t) wBc(t) uB vB wB              % velocity states
syms exd eyd ezd wxd wyd wzd uBd vBd wBd
syms q0c(t) q1c(t) q2c(t) q3c(t) q0 q1 q2 q3    % quaternions 
syms q0d q1d q2d q3d 
syms FBx FBy FBz TBx TBy TBz                    % summation of forces, torques
v1 = [xBc, diff(xBc, t), diff(xBc, t, t), yBc, diff(yBc, t), diff(yBc, t, t), ...
      zBc, diff(zBc, t), diff(zBc, t, t), exc, diff(exc, t), ...
      eyc, diff(eyc, t), ezc, diff(ezc, t), ...
      wxc, diff(wxc, t), wyc, diff(wyc, t), wzc, diff(wzc, t), ...
      q0c, diff(q0c, t), q1c, diff(q1c, t), q2c, diff(q2c, t), q3c, diff(q3c, t)];
v2 = [xB,  uB,           uBd,            yB,  vB,          vBd, ...
      zB,  wB,           wBd,            ex,  exd, ...
      ey,  eyd,          ez,  ezd, ...
      wx,  wxd,          wy,  wyd,       wz,  wzd, ...
      q0,  q0d,          q1,  q1d,       q2,  q2d,          q3,  q3d];

% rotation matrices - BCO style
% BCO = MakeBCO_euler(exc, eyc, ezc, 'XYZ');
% OCB = transpose(BCO);
% 
% make euler update equations - BCO
% OCBd = diff(OCB, t);
% wx_pe = simplify([0, 0, 1]*BCO*OCBd*[0; 1; 0]) - wx; % == 0
% wy_pe = simplify([1, 0, 0]*BCO*OCBd*[0; 0; 1]) - wy; % == 0
% wz_pe = simplify([0, 1, 0]*BCO*OCBd*[1; 0; 0]) - wz; % == 0
% wx_p1 = simplify(subs(wx_pe, v1, v2));
% wy_p1 = simplify(subs(wy_pe, v1, v2));
% wz_p1 = simplify(subs(wz_pe, v1, v2));
% [A1, b1] = equationsToMatrix([wx_p1, wy_p1, wz_p1], [exd, eyd, ezd]);
% eqnSet1 = simplify(A1\b1);  % exd = , eyd = , ezd = 

% rotation matrices - BCO quaternion
BCO = MakeBCO_quat(q0c, q1c, q2c, q3c);
OCB = transpose(BCO);

% make euler update equations - BCO
M = [  0, -wx, -wy, -wz; ...
      wx,   0,  wz, -wy; ...
      wy, -wz,   0,  wx; ...
      wz,  wy, -wz,   0];
eqnt = 0.5*M*[q0; q1; q2; q3] - [q0d; q1d; q2d; q3d];
eqnt0 = simplify(subs(eqnt(1), v1, v2));
eqnt1 = simplify(subs(eqnt(2), v1, v2));
eqnt2 = simplify(subs(eqnt(3), v1, v2));
eqnt3 = simplify(subs(eqnt(4), v1, v2));
[A1, b1] = equationsToMatrix([eqnt0, eqnt1, eqnt2, eqnt3], [q0d, q1d, q2d, q3d]);
eqnSet1 = simplify(A1\b1);  % q0d = , q1d = , q2d = , q3d =  

% positions
r_BO_B = MakeVec(xBc, yBc, zBc);    % Point B, generic on eraser
r_CMB_B = MakeVec(0, 0, 0);         % CM of eraser

% angular velocity
w_OB_B = MakeVec(wxc, wyc, wzc);

% linear velocities
O_v_BO_B = FirstTransport(r_BO_B, w_OB_B, t);
O_v_CMB_B = FirstTransport(r_CMB_B, w_OB_B, t);
O_v_CMO_B = O_v_CMB_B + O_v_BO_B;

% angular acceleration
a_OB_B = FirstTransport(w_OB_B, w_OB_B, t);

% linear accelerations
O_a_BO_B = SecondTransport(r_BO_B, w_OB_B, t);
O_a_CMB_B = SecondTransport(r_CMB_B, w_OB_B, t);
O_a_CMO_B = O_a_CMB_B + O_a_BO_B;

% inertia tensors
% Lx = length of side along x direction
% Ly = length of side along y direction
% Lz = length of side along z direction
I_B_B = (m/12)*MakeMat(Ly^2 + Lz^2, Lz^2 + Lx^2, Lx^2 + Ly^2, 0, 0, 0);

% angular momentum, rev0
% {OhOB,sys}B 
O_h_O_B_sys_B = I_B_B*w_OB_B;

% oddt OhOB
% {Od/dt OhOB,sys}B
Oddt_O_h_O_B_sys_B = FirstTransport(O_h_O_B_sys_B, w_OB_B, t) + CrossMe(O_v_BO_B, m*O_v_CMO_B);

% forces 
%   - define forces here
%   - requires forces defined here
Fg_O = [0; 0; -m*g];
Fg_B = BCO*Fg_O;
% %   - define forces in states file 
% %   - requires BCO and all there 
% sumF_B = MakeVec(FBx, FBy, FBz);

% torques
Tg_B = cross(r_CMB_B, Fg_B);    % gravity torque

% summation of forces
SumF_B = Fg_B - m*O_a_CMO_B;

% summation of torques
SumT_B = Tg_B - Oddt_O_h_O_B_sys_B;

% substitutions 
SumF_1_B = simplify(subs(SumF_B, v1, v2));
SumT_1_B = simplify(subs(SumT_B, v1, v2));

% convert to matrix
[A2, b2] = equationsToMatrix([SumF_1_B(1), SumF_1_B(2), SumF_1_B(3), ...
    SumT_1_B(1), SumT_1_B(2), SumT_1_B(3)], [uBd, vBd, wBd, wxd, wyd, wzd]);
eqnSet2 = simplify(A2\b2);    % eqns for [uBd, vBd, wBd, wxd, wyd, wzd]

% write to file for viewing
fid = fopen(eomsfile, 'w');
fprintf(fid, 'xd=, yd=, zd=, uBd=, vBd=, wBd=, q0d=, q1d=, q2d=, q3d=, wxd=, wyd=, wzd=\n');
% velocities
fprintf(fid, 'uB\n');
fprintf(fid, 'vB\n');
fprintf(fid, 'wB\n');
% accelerations 
fprintf(fid, '%s\n', eqnSet2(1));
fprintf(fid, '%s\n', eqnSet2(2));
fprintf(fid, '%s\n', eqnSet2(3));
% quat dots 
fprintf(fid, '%s\n', eqnSet1(1));
fprintf(fid, '%s\n', eqnSet1(2));
fprintf(fid, '%s\n', eqnSet1(3));
fprintf(fid, '%s\n', eqnSet1(4));
% rotation rate rates
fprintf(fid, '%s\n', eqnSet2(4));
fprintf(fid, '%s\n', eqnSet2(5));
fprintf(fid, '%s\n', eqnSet2(6));
fclose(fid);

% finish
disp('Done!');
    