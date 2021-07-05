function Generic_derive_and_simulation
% Function to simulate the motion of the eraser problem
% 
% TO DO
% 1. setup generic set for use with derivation, simulation, etc. 

clc
close all
clear
fclose('all');
clear('StatesFile');

% initial conditions 
% ICs = [x, y, z, u, v, w, q0, q1, q2, q3, p, q, r]
eomfile = 'generic_eoms.txt';
if exist(eomfile, 'file') ~= 2
    DeriveMe(eomfile);
end    

% test case 0 - no motion, no gravity
g = 0; m = 1; Lx = 1; Ly = 1; Lz = 1; ICs = [0, 0, 1.8, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]'; filename = 'generic_eoms_case0.csv'; tme = 0:0.1:10;

% parameters
p = [g, m, Lx, Ly, Lz];  % parameter vector




% simulate! :D - with events 
opts = odeset('Events', @(t, x)HitTheGroundStop(t, x, p));
[tme, outsB, te, ye, ie] = ode45(@(t, x)StatesFile(t, x, p, eomfile), tme, ICs, opts);

% % simulate! :D - without events
% [tme, outsB, te, ye, ie] = ode45(@(t, x)StatesFile(t, x, p, eomfile), tme, ICs);





% bring sim results to O frame
outsO = TranslateStates(tme, outsB);

% write time history to a file 
WriteToFile(filename, tme, outsB);

% plot things 
PlotMe(tme, p, outsB, outsO);

end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function [xd] = StatesFile(t, x, p, eomsname)
% Returns the state derivatives xd for a given time t, state vector x, and
% parameter set p.
% 
% INPUTS:
%   t = s, time
%   x = state vector
%   p = parameter vector = [g, m, Lx, Ly, Lz]
%   eomsname = file to read in for simulations
%
% OUTPUTS:
%   xd = state derivative vector

% unpack parameters
g = p(1);
m = p(2);
Lx = p(3);
Ly = p(4);
Lz = p(5);

% unpack states
% = [x, y, z, u, v, w, ex, ey, ez, p, q, r]
xB = x(1);
yB = x(2);
zB = x(3);
uB = x(4);
vB = x(5);
wB = x(6);
q0 = x(7);
q1 = x(8);
q2 = x(9);
q3 = x(10);
wx = x(11);
wy = x(12);
wz = x(13);

% read in equations for use
persistent eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8 eq9 eq10 eq11 eq12 eq13
if isempty(eq1)
    fid = fopen(eomsname, 'r');
    headline = fgetl(fid);   % header line
    eq1 = fgetl(fid);       % xd
    eq2 = fgetl(fid);       % yd
    eq3 = fgetl(fid);       % zd
    eq4 = fgetl(fid);       % xBdd
    eq5 = fgetl(fid);       % yBdd
    eq6 = fgetl(fid);       % zBdd
    eq7 = fgetl(fid);       % q0
    eq8 = fgetl(fid);       % q1
    eq9 = fgetl(fid);       % q2
    eq10 = fgetl(fid);      % q3
    eq11 = fgetl(fid);      % wxd
    eq12 = fgetl(fid);      % wyd
    eq13 = fgetl(fid);      % wzd
    fclose(fid);
end

% evaluate expressions
xd = NaN(length(x), 1);
xd(1) = eval(eq1);
xd(2) = eval(eq2);
xd(3) = eval(eq3);
xd(4) = eval(eq4);
xd(5) = eval(eq5);
xd(6) = eval(eq6);
xd(7) = eval(eq7);
xd(8) = eval(eq8);
xd(9) = eval(eq9);
xd(10) = eval(eq10);
xd(11) = eval(eq11);
xd(12) = eval(eq12);
xd(13) = eval(eq13);

end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function [outsO] = TranslateStates(tme, outsB)
% Converts states in B frame to O frame

% preallocate
outsO = NaN*outsB;

% rotate!
for i1 = 1:length(tme)
    
    % % euler angles
    % outsO(i1, 7:9) = outsB(i1, 7:9);
    % ex = outsB(i1, 7);
    % ey = outsB(i1, 8);
    % ez = outsB(i1, 9);
    % 
    % % rotation matrix
    % BCO = NaN(3, 3);
    % BCO(1, 1) = cos(ey)*cos(ez);
    % BCO(1, 2) = cos(ey)*sin(ez);
    % BCO(1, 3) = -sin(ey);
    % BCO(2, 1) = cos(ez)*sin(ex)*sin(ey) - cos(ex)*sin(ez);
    % BCO(2, 2) = cos(ex)*cos(ez) + sin(ex)*sin(ey)*sin(ez);
    % BCO(2, 3) = cos(ey)*sin(ex);
    % BCO(3, 1) = sin(ex)*sin(ez) + cos(ex)*cos(ez)*sin(ey);
    % BCO(3, 2) = cos(ex)*sin(ey)*sin(ez) - cos(ez)*sin(ex);
    % BCO(3, 3) = cos(ex)*cos(ey);
    % OCB = transpose(BCO);
    
    % quaternions
    outsO(i1, 7:10) = outsB(i1, 7:10);
    % rotation matrix 
    BCO = MakeBCO_quat(outsB(i1, 7), outsB(i1, 8), outsB(i1, 9), outsB(i1, 10));
    OCB = transpose(BCO);
    
    % positions
    pos_B = [outsB(i1, 1); outsB(i1, 2); outsB(i1, 3)];
    outsO(i1, 1:3) = transpose(OCB*pos_B);
    
    % velocities
    vel_B = [outsB(i1, 4); outsB(i1, 5); outsB(i1, 6)];
    outsO(i1, 4:6) = transpose(OCB*vel_B);
    
    % angular velocities
    ang_B = [outsB(i1, 10); outsB(i1, 11); outsB(i1, 12)];
    outsO(i1, 10:12) = transpose(OCB*ang_B);
    
end

end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function PlotMe(tme, p, outsB, outsO)

% plot the results - O frame
mrkrs = {'-o'; '-^'; '-s'};
figure('Color', 'w');
subplot(4, 1, 1);
hold on
plot(tme, outsO(:, 1), mrkrs{1});
plot(tme, outsO(:, 2), mrkrs{2});
plot(tme, outsO(:, 3), mrkrs{3});
ylabel('Position [m]', 'interpreter', 'latex');
grid on
legend({'$\vec{i}_{\bar{O}}$', '$\vec{j}_{\bar{O}}$', '$\vec{k}_{\bar{O}}$'}, 'location', 'southwest', 'interpreter', 'latex');

subplot(4, 1, 2);
hold on
plot(tme, outsO(:, 4), mrkrs{1});
plot(tme, outsO(:, 5), mrkrs{2});
plot(tme, outsO(:, 6), mrkrs{3});
ylabel('Velocity [m]', 'interpreter', 'latex');
grid on

subplot(4, 1, 3);
hold on
plot(tme, outsO(:, 7)*(180/pi), mrkrs{1});
plot(tme, outsO(:, 8)*(180/pi), mrkrs{2});
plot(tme, outsO(:, 9)*(180/pi), mrkrs{3});
ylabel('Euler angles [deg]', 'interpreter', 'latex');
grid on

subplot(4, 1, 4);
hold on
plot(tme, outsO(:, 10), mrkrs{1});
plot(tme, outsO(:, 11), mrkrs{2});
plot(tme, outsO(:, 12), mrkrs{3});
ylabel('Angular rates [deg/s]', 'interpreter', 'latex');
grid on
xlabel('Time [s]', 'interpreter', 'latex');

% plot the results - B frame
mrkrs = {'-o'; '-^'; '-s'};
figure('Color', 'w');
subplot(4, 1, 1);
hold on
plot(tme, outsB(:, 1), mrkrs{1});
plot(tme, outsB(:, 2), mrkrs{2});
plot(tme, outsB(:, 3), mrkrs{3});
ylabel('Position [m]', 'interpreter', 'latex');
grid on
legend({'$\vec{i}_{\bar{B}}$', '$\vec{j}_{\bar{B}}$', '$\vec{k}_{\bar{B}}$'}, 'location', 'southwest', 'interpreter', 'latex');

subplot(4, 1, 2);
hold on
plot(tme, outsB(:, 4), mrkrs{1});
plot(tme, outsB(:, 5), mrkrs{2});
plot(tme, outsB(:, 6), mrkrs{3});
ylabel('Velocity [m]', 'interpreter', 'latex');
grid on

subplot(4, 1, 3);
hold on
plot(tme, outsB(:, 7)*(180/pi), mrkrs{1});
plot(tme, outsB(:, 8)*(180/pi), mrkrs{2});
plot(tme, outsB(:, 9)*(180/pi), mrkrs{3});
ylabel('Euler angles [deg]', 'interpreter', 'latex');
grid on

subplot(4, 1, 4);
hold on
plot(tme, outsB(:, 10), mrkrs{1});
plot(tme, outsB(:, 11), mrkrs{2});
plot(tme, outsB(:, 12), mrkrs{3});
ylabel('Angular rates [deg/s]', 'interpreter', 'latex');
grid on
xlabel('Time [s]', 'interpreter', 'latex');

end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function WriteToFile(filename, tme, outsB)

fid = fopen(filename, 'w');
if length(outsB(1, :)) == 12
    fprintf(fid, '     t(s),    xB(m),    yB(m),    zB(m),  uB(m/s),  vB(m/s),  wB(m/s),  ex(rad),  ey(rad),  ez(rad), p(rad/s), q(rad/s), r(rad/s)\n');
else
    fprintf(fid, '     t(s),    xB(m),    yB(m),    zB(m),  uB(m/s),  vB(m/s),  wB(m/s),       q0,       q1,       q2,       q3, p(rad/s), q(rad/s), r(rad/s)\n');
end
for i1 = 1:length(tme)
    strg = sprintf('%9.4f', tme(i1));
    for j1 = 1:length(outsB(1, :))
        strg = sprintf('%s,%9.4f', strg, outsB(i1, j1));
    end
    fprintf(fid, '%s\n', strg);
end
fclose(fid);

end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function [val, ister, dire] = HitTheGroundStop(t, x, p)
% % rotation matrix 
% BCO = MakeBCO(x(7), x(8), x(9), 'XYZ');
% OCB = transpose(BCO);

% rotation matrix 
BCO = MakeBCO_quat(x(7), x(8), x(9), x(10));
OCB = transpose(BCO);

% calculate ground error
rBO_O = OCB*[x(1); x(2); x(3)];
val = rBO_O(3);
ister = 1;
dire = 0;
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function DeriveMe(eomsfile)
% Generic derivation outline 
% 
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

end