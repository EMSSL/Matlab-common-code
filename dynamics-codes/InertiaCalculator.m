%{
code to calculate moments of inertia about n rigid bodies
christopher d. yoder
10/27/2021

How it works:
1. Specify each rigid body in a matrix ahead of time. First row is first
rigid body, second row(body), etc. xyz body axes must align with the XYZ 
coordinate frame used to specify locations of xCM, yCM, zCM. Columns are:
    m, xCM, yCM, zCM, Ixx, Iyy, Izz, Ixy, Ixz, Iyz
2. Everything is defined with respect to XYZ=(0,0,0) in global space. 
3. Ixx, etc must be defined about the body CM, specified as xCM, yCM, zCM. 
%}

% preamble
clc
close all
clear 
fclose('all');

% bodies matrix: T-handle problem
% Rod 1: slender along iO, kO is up, jO is RHR. Length L1, radius R1. CM at (0,0,0). 
% Rod 2: slender along jO, kO is up, iO is RHR. Length L2, radius R2. CM at (0,L2/2,0). 
% m1 = 0.1; L1 = 0.1; r1 = 0.01; m2 = m1; L2 = L1; r2 = r1;
syms m1 m2 L1 L2 r1 r2
bodymatrix = [m1, 0, 0, 0, m1*r1*r1/2, m1*(3*r1*r1 + L1*L1)/12, m1*(3*r1*r1 + L1*L1)/12, 0, 0, 0;...
    m2, 0, L2/2, 0, m2*(3*r2*r2 + L2*L2)/12, m2*r2*r2/2, m2*(3*r2*r2 + L2*L2)/12, 0, 0, 0];

% % bodies matrix: point mass sanity check
% % Point 1: CM at (1,0,0), m1 = 1
% % Point 2: CM at (1,1,0), m2 = 1
% % Point 3: CM at (1,1,1), m3 = 1
% % Point 4: CM at (1,0,1), m4 = 1
% bodymatrix = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0; ...
%     1, 1, 1, 0, 0, 0, 0, 0, 0, 0; ...
%     1, 1, 1, 1, 0, 0, 0, 0, 0, 0; ...
%     1, 1, 0, 1, 0, 0, 0, 0, 0, 0];





% calculate CM position
xCM_O = sum(bodymatrix(:, 1).*bodymatrix(:, 2))/sum(bodymatrix(:, 1));
yCM_O = sum(bodymatrix(:, 1).*bodymatrix(:, 3))/sum(bodymatrix(:, 1));
zCM_O = sum(bodymatrix(:, 1).*bodymatrix(:, 4))/sum(bodymatrix(:, 1));
rCMO_O = [xCM_O; yCM_O; zCM_O];
disp('CM position:');
disp(rCMO_O);
if isnumeric(xCM_O) == 1
    fprintf("xCM_O = %0.4f;\n", xCM_O);
    fprintf("yCM_O = %0.4f;\n", yCM_O);
    fprintf("zCM_O = %0.4f;\n\n", zCM_O);
else
    fprintf("xCM_O = %s;\n", xCM_O);
    fprintf('yCM_O = %s;\n', yCM_O);
    fprintf("zCM_O = %s;\n\n", zCM_O);
end


% calculate inertia tensor 
ICM = [0, 0, 0; 0, 0, 0; 0, 0, 0];
for i1 = 1:length(bodymatrix(:, 1))
    % get distance vector
    d = bodymatrix(i1, 2:4).' - rCMO_O;
    % make inertia tensor about CM_i
    Ixx = bodymatrix(i1, 5);
    Iyy = bodymatrix(i1, 6);
    Izz = bodymatrix(i1, 7);
    Ixy = bodymatrix(i1, 8);
    Ixz = bodymatrix(i1, 9);
    Iyz = bodymatrix(i1, 10);
    ICM_i = [Ixx, Ixy, Ixz; Ixy, Iyy, Iyz; Ixz, Iyz, Izz];
    % move to big CM
    ICM = ICM + ICM_i + d*(d.');
end
disp('ICM tensor:');
disp(ICM);
if isnumeric(ICM(1,1)) == 1
    fprintf("Ixx = %0.4f;\n", ICM(1,1));
    fprintf("Iyy = %0.4f;\n", ICM(2,2));
    fprintf("Izz = %0.4f;\n", ICM(3,3));
    fprintf("Ixy = %0.4f;\n", ICM(1,2));
    fprintf("Ixz = %0.4f;\n", ICM(1,3));
    fprintf("Iyz = %0.4f;\n", ICM(2,3));
else
    fprintf("Ixx = %s;\n", ICM(1,1));
    fprintf("Iyy = %s;\n", ICM(2,2));
    fprintf("Izz = %s;\n", ICM(3,3));
    fprintf("Ixy = %s;\n", ICM(1,2));
    fprintf("Ixz = %s;\n", ICM(1,3));
    fprintf("Iyz = %s;\n", ICM(2,3));
end