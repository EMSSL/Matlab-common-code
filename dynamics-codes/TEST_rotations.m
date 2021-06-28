% testing rotations, etc. 
% christopher d. yoder
% 06/27/2021

clc
close all
clear 
fclose('all');
% UsePackage('remove', 'DynPackage');
r2d = 180/pi;
d2r = 1/r2d;

% test 1: derive components from BCO matrix
% -------------------------------------------------
seq = 'XYZ';  
a1 = rand()*2*pi - pi;
a2 = rand()*2*pi - pi;
a3 = rand()*2*pi - pi;
BCO = MakeBCO_euler(a1, a2, a3, seq)                            % make BCO "truth"
[t1, t2, t3, q0, q1, q2, q3] = MatrixToComponents(BCO, seq);    % make components
BCO_e = MakeBCO_euler(t1, t2, t3, seq);
BCO_q = MakeBCO_quat(q0, q1, q2, q3);

% compare accuracies
BCO - BCO_e         % difference in truth and estimate using Eulers
norm(BCO - BCO_e)   % norm in truth and estimate using Eulers
BCO - BCO_q         % difference in truth and estimate using Quats
norm(BCO - BCO_q)   % norm in truth and estimate using Eulers