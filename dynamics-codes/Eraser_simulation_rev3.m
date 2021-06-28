function Eraser_simulation_rev3
% Function to simulate the motion of the eraser problem

clc
close all
clear
fclose('all');
clear('StatesFile');

% initial conditions 
% ICs = [x, y, z, u, v, w, q0, q1, q2, q3, p, q, r]
eomfile = 'eraser_equations_rev3.txt';



% test case 0 - no motion, no gravity
g = 0; m = 1; Lx = 1; Ly = 1; Lz = 1; ICs = [0, 0, 1.8, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]'; filename = 'case0q.csv'; tme = 0:0.1:10;

% test case 1 - translation
% g = 1.625; m = 1; Lx = 1; Ly = 1; Lz = 1; ICs = [0, 0, 1.8, 10, 0, 10, 1, 0, 0, 0, 0, 0, 0]'; filename = 'case1q.csv'; tme = 0:0.1:10;

% test case 2 - spin only
% g = 1.625; m = 1; Lx = 1; Ly = 1; Lz = 1; ICs = [0, 0, 1.8, 0, 0, 0, 1, 0, 0, 0, 0, 3.14, 0]'; filename = 'case2q.csv'; tme = 0:0.1:10;

% test case 3 - altitude
% g = 1.625; m = 1; Lx = 1; Ly = 1; Lz = 1; ICs = [0, 0, 1.8, 0, 0, 10, 1, 0, 0, 0, 0, 0, 0]'; filename = 'case3q.csv'; tme = 0:0.1:10;

% test case 4 - curve ball
% g = 1.625; m = 1; Lx = 1; Ly = 1; Lz = 1; ICs = [0, 0, 1.8, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0]'; filename = 'case4q.csv'; tme = 0:0.1:10;

% test case 5 - drop test
% g = 1.625; m = 1; Lx = 1; Ly = 1; Lz = 1; ICs = [0, 0, 1.8, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0]'; filename = 'case5q.csv'; tme = 0:0.1:10;


% simulate! :D
p = [g, m, Lx, Ly, Lz];  % parameter vector
opts = odeset('Events', @(t, x)HitTheGroundStop(t, x, p));
[tme, outsB, te, ye, ie] = ode45(@(t, x)StatesFile(t, x, p, eomfile), tme, ICs, opts);

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
fprintf(fid, '     t(s),    xB(m),    yB(m),    zB(m),  uB(m/s),  vB(m/s),  wB(m/s),  ex(rad),  ey(rad),  ez(rad), p(rad/s), q(rad/s), r(rad/s)\n');
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