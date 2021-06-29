% code to recreate EDO hw using PSO 

clc
close all
clear

% problem 
popSize = 25;               % population size
iterMin = 10;               % minimum iteration count before kickout
iterMax = 1e4;              % maximum iterations
reqdIter = 50;              % required iterations for convergence
thetaMax = 0.83;            % max theta value
thetaMin = 0.4;             % min theta value
cLocal = 2.4;               % local weight factor
cGlobal = 1.7;              % global weight factor
cVec = [cLocal, cGlobal];
moveLims = [3, 15];

% other constants and whatnot
x0 = [10, 20, 30, 40, 50, 0, 0, 0, 0, 0];     % initial values

% optimize!
objFun = 'hw3p5objF';
[xStar, fStar, iter] = ParticleSwarm(objFun, x0, 'display', 'exit', ...
    'iterMax', 1e4, 'theta', 0, 'moveLim', 10);