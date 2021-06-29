% ParticleSwarm(objFun, x0) uses the particle swarm optimization formulation 
% to estimate the optimuym solution to the problem in objFun using the
% design variables in x0 as a starting location. 
% [xStar, fStar, iter] = ParticleSwarm(objFun, x0) returns the optimal
% solution xStar, the function value fStar, and the iteration count to
% generate the solution. 
%   -----------------------------------------------------------------
% INPUTS:
%   objFun    =    function handle with the objective function to optimize.
%                  Call the function name with the '@objFun' convention
%                  from outside the PSO file to allow extra parameters to
%                  be passed to functions.
%   x0        =    initial design variable string
%   -----------------------------------------------------------------
% OUTPUTS:
%   xStar     =    optimal design variables for the problem
%   fStar     =    optimal function value for the problem
%   iter      =    iteration number to obtain the solution
%   feval     =    number of function evaluations
%   -----------------------------------------------------------------
% OPTIONS: {varargin}
%   PROBLEM FORMULATION:
%   LB        =    lower bounds for x vector used for initial populations
%   UB        =    upper bounds for x vector used for initial populations
%   -----------------------------------------------------------------
%   SOLVER PARAMETERS:
%   display   =    (default is 'off') displays information. Options are:
%                   'exit' to display exit message and final solution
%                   'iter' to display the results of each iteration - no go
%   popSize   =    (default is 25) number of designs used to optimize.
%   reqdIter  =    (default is 50) number of iterations for which the
%                  global optimal solution must be constant for
%                  convergence.
%   iterMin   =    (default is 10) minimum number of iterations before
%                   algorithm can return best solution
%   iterMax   =    (default is 1e5) maximum number of iterations before
%                   algorithm must return best solution, regardless of
%                   optimality
%   theta     =    (default is 'dynamic', range of [0.4, 0.9]). Theta can
%                   be set to a static value by providing a scalar. The dynamic range can
%                   be changed by providing a vector. 
%   cVec      =    (default is [2.4, 1.7]) values of [cLocal, cGlobal] 
%                   used for the calculations of the velocities.
%   moveLim   =    (default is +/- 1.5) move restriction for x values
%   waitbar   =    (default is 'off'|'on') used for objective function 
%                   runs which are long. 
%   golden    =     (default is 'off'|'on') uses goldenSection algorithm to
%                   converge quickly to a solution. Trigger criteria (10) is
%                   when one solution continues to be the best while making
%                   small improvements. 
%   trigger   =     (default is 10) 
%   -----------------------------------------------------------------
% NOTES:
% 1. The objective function for PSO generally is either unconstrained and
%       unbounded OR contains a built-in penalty function to correct for
%       infeasible solutions. 

% Written Fall 2016 for EDO course
% Modifications:
% 1. CDY, 6/29/17, generalized the function to be used by any function
% 2. CDY, 7/7/18,  fixed objective function problems, fixed LB and UB
%                  problems in the code, added feval counter, updated the
%                  display 
%
% To Do:
% 3. CDY, 7/26/19, Added in feval to the output list

function [xStar, fStar, iter, feval] = ParticleSwarm(objFun, x0, varargin)

% STEP 1: Define the defaults input output stuffs
popSize = 25;           % population size
iterMin = 1;           % minimum iteration count before kickout
iterMax = 1e5;          % maximum iterations
reqdIter = 50;          % required iterations for convergence
thetaFlag = 0;
thetaMax = 0.9;         % max theta value
thetaMin = 0.4;         % min theta value
cLocal = 2.4;           % local weight factor
cGlobal = 1.7;          % global weight factor
moveLim = 1.5;          % move limit x
displayFlag = 0;        % 0 = 'off', 1 = 'exit', 2 = 'iter', 3 = '2Dplot', 4 = '3Dplot'
bestTOL = 1e-6;
lbFlag = 0;
ubFlag = 0;
LB = [];
UB = [];
feval = 0;
waitbarFlag = 0;
posbar = [];
goldenCntr = 0;
goldenTrigger = 10;     
goldenFlag = 0;
goldenIndx = 0;
goldenOld = goldenIndx;

% override default values if user specifies it
if isempty(varargin) ~= 1
    m1 = length(varargin)/2;
    for i = 1:m1
        switch char(varargin(2*(i - 1) + 1))
            case 'popSize'
                popSize = cell2mat(varargin(2*(i - 1) + 2));    % change popSize
                
            case 'iterMin'  
                iterMin = cell2mat(varargin(2*(i - 1) + 2));    % change iterMin
                
            case 'iterMax'
                iterMax = cell2mat(varargin(2*(i - 1) + 2));    % change iterMax
                
            case 'reqdIter'
                reqdIter = cell2mat(varargin(2*(i - 1) + 2));   % change reqdIter
                
            case 'theta'
                
                % set thetas
                thetaIn = cell2mat(varargin(2*(i - 1) + 2));    % get theta input
                if length(thetaIn) == 1          % if theta is scalar
                    thetaFlag = 1;               % set flag 1
                    
                elseif length(thetaIn) == 2      % if theta is vector
                    thetaFlag = 0;               % set flag 0
                    
                    if thetaIn(2) < thetaIn(1)   % error handling
                        error('theta(2) must be greater than theta(1).');
                    end
                    
                    thetaMin = thetaIn(1);       % set theta min
                    thetaMax = thetaIn(2);       % set theta max
                    
                end

            case 'cVec'
                cVec = cell2mat(varargin(2*(i - 1) + 2));    % get c vector input
                cLocal = cVec(1);
                cGlobal = cVec(2);
                
            case 'moveLim'
                moveLim = cell2mat(varargin(2*(i - 1) + 2)); % move limits for initial populations
                
            case 'LB'
                LB = cell2mat(varargin(2*(i - 1) + 2)); % move limits for initial populations
                lbFlag = 1;
                
            case 'UB'
                UB = cell2mat(varargin(2*(i - 1) + 2)); % move limits for initial populations
                ubFlag = 1;
                
            case 'display'
                
                % set new flags accordingly
                switch char(varargin(2*(i - 1) + 2))
                    case 'off'
                        displayFlag = 0;
                    case 'exit'
                        displayFlag = 1;
                        disp('Beginning Particle Swarm Optimization method. . .');
                    case 'iter'
                        displayFlag = 2;
                        disp('Beginning Particle Swarm Optimization method. . .');
                end
                
            case 'waitbar'
                if strcmp(varargin(2*(i - 1) + 2), 'on') == 1
                    waitbarFlag = 1;
                end
                
            case 'golden'
                switch char(varargin(2*(i - 1) + 2))
                    case 'off'
                        goldenFlag = 0;
                    case 'on'
                        goldenFlag = 1;
                    otherwise
                        error('golden argument not recognized.');
                end
            case 'trigger'
                goldenTrigger = varargin(2*(i - 1) + 2);
               
        end
    end
end

% handle bounding errors
if lbFlag + ubFlag > 0 && lbFlag + ubFlag ~= 2
    error('Both LB and UB must be specified.');
end
if length(LB) ~= length(UB)
    error('Both LB and UB must have the same number of elements.');
end
if isrow(LB) ~= 1
    LB = LB';       % for math later
end
if isrow(UB) ~= 1
    UB = UB';       % for math later
end

% handle inputs
if isrow(x0) ~= 1
    x0 = transpose(x0);
end


% STEP 2: Initialize population with temperature
m2 = length(x0);                                % length of x vector
popString = NaN(popSize, 1 + m2);               % preallocate
randValsX = NaN(popSize, m2);                   % preallocate
initVals = x0;                                  % remember initial values
if lbFlag + ubFlag == 2
    for k3 = 1:popSize
        randValsX(k3, :) = (UB - LB).*rand(1, m2) + LB;
    end
else
    randValsX = (rand(popSize, m2) - 0.5)*moveLim;  % generate random populations
end

% SCHEDULE - popString(ith, :) = [x2, x3, x4, x5, y2, y3, y4, y5, Feval]
if waitbarFlag == 1
    wb = waitbar(0, 'Running objective function...');
    posbar = get(wb, 'Position');
end
for k1 = 1:popSize                                                          % for all members in population
    popString(k1, 1:m2) = initVals + randValsX(k1, :);                      % generate initial positions
    popString(k1, 1:m2) = BoundPopulation(popString(k1, 1:m2), LB, UB);     % adjust for bounds
    
    % popString(k1, 11) = eval([objFun, '(transpose(popString(k1, 1:m2)))']); % find function value  
    popString(k1, m2 + 1) = objFun(popString(k1, 1:m2));    
    feval = feval + 1;
    if waitbarFlag == 1
        waitbar(k1/popSize, wb, 'Running objective function...');
    end
    
end
if waitbarFlag == 1
    posbar = get(wb, 'Position');
    close(wb);
end

% set velocities 
cntr = 1;                                           % set iteration counter to 1
convergeCounter = 0;                                % set convergence counter to 0
velVals = zeros(popSize, m2);                       % set velocity values to zero
[~, minFunInd] = min(popString(:, m2 + 1));             % find global best init fun value
globalBest = popString(minFunInd, :);               % populate the design string
localBest = popString;                              % initialize local best matrix

% display information headers if required
if displayFlag == 2
    % fprintf('%-10s\t', '%-10s\t', '%-10s\t',  '%-10s\n', 'Iter', 'fVal*', 'ConvCntr', 'Best#');
    fprintf('%-9s\t %-9s\t %-9s\t %-9s\t %-9s\t %-9s\n', 'Iter', 'Global', 'dfVal', 'Converge', 'Best', 'Function');
    fprintf('%-9s\t %-9s\t %-9s\t %-9s\t %-9s\t %-9s\n', ' ', 'Best', 'diter ', 'Counter', 'Particle', 'Calls');
    fprintf('---------    ---------   ---------   ---------   ---------   ---------  \n');
end


% iterate!
exitFlag = 0;
while exitFlag == 0
    
    % initialize waitbar
    if waitbarFlag == 1
        wb = waitbar(0, 'Running objective function...', 'Position', posbar);
    end

    % STEP 3: get best and calculate velocity
    if thetaFlag == 1               % if theta is scalar
        theta = thetaIn;            % rename and move on
    else                            % if vector, calculate theta using line
        theta = thetaMax - cntr*(thetaMax - thetaMin)/iterMax;
    end
    
    % caluclate new velocities
    for k2 = 1:popSize                                                              % for all members of population
        velVals(k2, :) = theta*velVals(k2, :) + ...                                 % calculate new velocity
            cLocal*rand(1, 1)*(localBest(k2, 1:m2) - popString(k2, 1:m2)) + ...
            cGlobal*rand(1, 1)*(globalBest(1, 1:m2) - popString(k2, 1:m2));
        popString(k2, 1:m2) = popString(k2, 1:m2) + velVals(k2, :);                 % generate new population string
        
%         % CHECK TO MAKE SURE THAT THE DESIGN STRING IS WITHIN THE BOUNDS AS
%         % STATED BY LB AND UB
%         keyboard
%         if lbFlag + ubFlag == 2
%             if sum(popString(k2, 1:m2) > UB) > 0                % if any values are above the UB
%                 badUp = find(popString(k2, 1:m2) > UB);         % find which elements they are
%                 for k4 = 1:length(badUp)                        % for all bad elements
%                     popString(k2, badUp(k4)) = UB(badUp(k4));   % replace the element with the upper bound
%                 end
%             end
%             if sum(popString(k2, 1:m2) < LB) > 0                % if any values are below the LB
%                 badLo = find(popString(k2, 1:m2) < LB);         % find which elements they are
%                 for k4 = 1:length(badLo)                        % for all bad elements
%                     popString(k2, badLo(k4)) = UB(badLo(k4));   % replace the element with the lower bound
%                 end
%             end
%         end

        popString(k2, 1:m2) = BoundPopulation(popString(k2, 1:m2), LB, UB);
        % popString(k2, m2 + 1) = eval([objFun, '(transpose(popString(k2, 1:m2)))']); % find function value
        popString(k1, m2 + 1) = objFun(popString(k2, 1:m2));  
        feval = feval + 1;
        if waitbarFlag == 1
            waitbar(k2/popSize, wb, 'Running objective function...');
        end
    end
    if waitbarFlag == 1
        posbar = get(wb, 'Position');
        close(wb);
    end


    
    % STEP 4: update the local and global bests
    convergeIncrFlag = 1;
    bestParticle = NaN;
    dfVal = NaN;
    
    for k3 = 1:popSize

        % update local bests
        if popString(k3, m2 + 1) < localBest(k3, m2 + 1)
            localBest(k3, :) = popString(k3, :);
        end

        % update global best
        if popString(k3, m2 + 1) < globalBest(m2 + 1) - bestTOL 
            dfVal = abs(globalBest(m2 + 1)) - abs(popString(k3, m2 + 1));
            globalBest = popString(k3, :);
            convergeIncrFlag = 0;
            bestParticle = k3;
            goldenIndx = k3;
            if goldenFlag == 1
                if goldenIndx == goldenOld
                    goldenCntr = goldenCntr + 1;
                else
                    goldenOld = goldenIndx;
                    goldenCntr = 1;
                end
            end
        end

    end
    
    if displayFlag == 2 
%         if isnan(bestParticle) ~= 1 || cntr == 1
%             fprintf('%-9d\t %-9.4e\t %-9.4e\t %-9d\t %-9d\t %-9.4e\n', ...
%                 cntr, globalBest(end), dfVal, convergeCounter, bestParticle, feval);
%         end
        fprintf('%-9d\t %-9.4e\t %-9.4e\t %-9d\t %-9d\t %-9.4e\n', ...
            cntr, globalBest(end), dfVal, convergeCounter, bestParticle, feval);
    end

    % count for the same stuffs over and over
    convergeCounter = convergeIncrFlag*(convergeCounter + convergeIncrFlag);
    cntr = cntr + 1;
    
    % monitor exit flag
    if convergeCounter > reqdIter 
        exitFlag = 1;
        msge = ['PSO exited normally due to convergence counter being greater', ...
                ' than the required iteration count and less than iterMax.'];
    end
    if cntr > iterMax
        exitFlag = 2;
        msge = ['PSO exited because the iteration counter exceeded the ', ...
            ' allowable iteration count, iterMax.'];
    end 
    if goldenCntr > goldenTrigger
        exitFlag = 3;
        msge = ['PSO exited because the golden section counter exceeded the ', ...
            ' golden iteration count, goldenTrigger.'];
    end



end 

% run golden section
if goldenFlag == 1 && exitFlag == 3
    % run powells method
    [xStar, fStar, ~, ~] = PowellsMethod(objFun, globalBest(1:m2));
else
    % report outputs
    xStar = globalBest(1:m2);
    fStar = globalBest(m2 + 1);
    iter = cntr;
end

% display exit results
if displayFlag == 1 || displayFlag == 2
    disp(' ');
    disp(msge);
    disp('Optimal design values: ');
    disp(xStar);
    disp('Optimal design function value:');
    disp(fStar);
    disp('Iteration count:');
    disp(cntr);
    disp('Function calls:');
    disp(feval);
end

end




function [PopStringNew] = BoundPopulation(PopStringOld, LB, UB)
    
% CHECK TO MAKE SURE THAT THE DESIGN STRING IS WITHIN THE BOUNDS AS
% STATED BY LB AND UB
PopStringNew = PopStringOld;                            % preallocate as old string
if isempty(UB) ~= 1
    if sum(PopStringOld > UB) > 0                       % if any values are above the UB
        badUp = find(PopStringOld > UB);                % find which elements they are
        for k4 = 1:length(badUp)                        % for all bad elements
            PopStringNew(badUp(k4)) = UB(badUp(k4));    % replace the element with the upper bound
        end
    end
end
if isempty(LB) ~= 1
    if sum(PopStringOld < LB) > 0                       % if any values are below the LB
        badLo = find(PopStringOld < LB);                % find which elements they are
        for k4 = 1:length(badLo)                        % for all bad elements
            PopStringNew(badLo(k4)) = LB(badLo(k4));    % replace the element with the lower bound
        end
    end
end

end