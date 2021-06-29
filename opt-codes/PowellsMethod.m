function [xopt, Fopt, iter, funEvals] = PowellsMethod(objFun, x0, varargin)
% [xopt, Fopt, iter] = PowellsMethod(objFun, x0, varargin) returns the
% optimal design string, function evaluation, and iteration count for the
% objFun function with starting point x0 using Powell's method
%   -----------------------------------------------------------------
% INPUTS:
%   objFun  = Function handle to the objective function to be optimized.
%             objFun should return a scalar output. Constraints and bounds
%             should be handled using a penalty function and
%             pseudoobjective function method within objFun itself. 
%   x0      = Starting point for Powell's method. 
%   -----------------------------------------------------------------
% OUTPUTS:
%   xopt        = Optimum design string.
%   Fopt        = Function value at xopt.
%   iter        = Iteration count 
%   funEvals    = number of function evaluations.
%   -----------------------------------------------------------------
% OPTIONS:
%   'display'   =   'off' (default) Display the information during the
%                   progress of the algorithm. Can be set to 'iter' to see
%                   iteration information or 'soln' to see results upon
%                   exit. 
%   'LB'        =   -INF (default). Lower bounds on design string elements
%   'UB'        =   +INF (default). Upper bounds on design string elements
%   'domain'    =   'continuous' (default). For discrete problems, set
%                   'domain' to 'discrete'. The parameter 'incr' must also
%                   be set. 
% CONTINUOUS OPTIONS:
%   'delta'     =   0.01 (default) Step size for bounding phase algorithm. 
%   'tol'       =   1e-5 (default) Solution tolerance to exit the
%                   algorithm. Measured as the differene in objective 
%                   function value between iterations 
%   'iterMax'   =   1e2 (default) Maximum iterations before returning the
%                   current solution. 
% DISCRETE OPTIONS:
%   'incr'      =   1 (default). 'incr' is the
%                   increment by which a discrete design variable is 
%                   increased or decreased.
%       Example: Design string is: [car engine, gear box, fuel tank]
%                Choices are: five engines, two gear boxes, and three tanks
%                Set 'incr' to [1, 1, 1], set 'LB' to [1, 1, 1], set 'UB'
%                to [5, 2, 3]. Then, the design string [3, 2, 2] would
%                represent [Car engine 3, Gear box 2, Fuel tank 2].
%       Example: Design string is: [diameter, thread length, color]
%                Bolt diameter could be [0.25, 0.5, 0.75, 1]
%                Thread length could be [0.5, 1, 1.5, 2]
%                Color could be ['red', 'blue', 'yellow']
%                Set 'incr' to [0.25, 0.5, 1], set 'LB' to [0.25, 0.5, 1], 
%                set 'UB' to [1, 2, 3]. Then, the design string [3, 2, 2]
%                would represent [0.75, 1, 'blue'].
%                

% Written by Christopher D. Yoder in Fall 2016
% Modifications:
% 1. CDY, 6/30/17, Updated input syntax to rearrange to match other
%    functions, added varargin capabilities, updated help.
% 2. CDY, 7/6/2017, Added bounds to design string to pass to boundingPhase
%  

% To Do
% 1. Add in several stopping critera 
%   - max function evals
%   - change in function value tolerance
% 2. Add bounding to the algorithm


% STEP 1: Initialize stuffs
sMat = eye(length(x0));
tolerance = 1e-5;
iter = 0;
delta = 0.01;
exitFlag = 0;
dispFlag = 0;
iterMax = 1e2;
funEvals = 1;
dFunTol = 1e-2;
LB = -inf*ones(1, length(x0));
UB = inf*ones(1, length(x0));

% input handling
if isempty(varargin) ~= 1
    for i1 = 1:length(varargin)/2
        switch varargin{2*(i1 - 1) + 1}
            case 'tol'
                tolerance = varargin{2*(i1 - 1) + 2};
            case 'display'
                
                switch varargin{2*(i1 - 1) + 2}
                    case 'iter'
                        dispFlag = 1;
                    case 'soln'
                        dispFlag = 2;
                    case 'off'
                        dispFlag = 0;
                    otherwise
                        error('Display flag not recognized.');
                end
            case 'delta'
                delta = varargin{2*(i1 - 1) + 2};
            case 'iterMax'
                iterMax = varargin{2*(i1 - 1) + 2};
            case 'LB'
                LB = varargin{2*(i1 - 1) + 2};
                if length(LB) ~= length(x0)
                    error(['Length of LB (', num2str(length(LB)), ...
                        ') does not equal Length of x0 (', ...
                        num2str(length(x0))]);
                end
            case 'UB'
                UB = varargin{2*(i1 - 1) + 2};
                if length(UB) ~= length(x0)
                    error(['Length of UB (', num2str(length(UB)), ...
                        ') does not equal Length of x0 (', ...
                        num2str(length(x0))]);
                end
            otherwise
                error(['Option "', varargin{2*(i1 - 1) + 1}, '" is not recognized.']);
        end
    end
end
                

% error handling
if isrow(x0) == 1
    x0 = x0';
end

% disp initial
if dispFlag == 1
    fprintf('%-10s\t %-15s\t %-15s\t %-15s\t \n', 'Iteration', 'Function Count', 'Function Value', 'A*');
    fprintf('______________________________________________________________ \n');
end
                    

% STEP 2: Powells method: find alpha star 
fopt1 = objFun(x0);
while exitFlag == 0
    
    x_orig = x0;
    iter = iter + 1;
    for i1 = 1:length(x0)

        % bounding phase
        s0 = sMat(:, i1);
        % opts = 'iter';
        % opts = 'none';

%         [xLEFT, xRIGHT, ~, ~, ~] = boundingPhase(x0, delta, ...
%             objFun, 'direction', s0, 'display', opts);
%         [xLEFT, xRIGHT, ~, ~, ~, funs] = boundingPhase(objFun, x0, delta, ...
%             'direction', s0, 'display', 'none');
        [xLEFT, xRIGHT, ~, ~, ~, funs] = boundingPhase(objFun, x0, delta, ...
            'direction', s0, 'display', 'none', 'LB', LB, 'UB', UB);
        funEvals = funEvals + funs;

        % golden section phase for alpha star
        % [xopt, ~, ~] = goldenSection(xLEFT, xRIGHT, objFun);
        [xopt, ~, ~, funs] = goldenSection(objFun, xLEFT, xRIGHT);
        x0 = xopt';
        funEvals = funEvals + funs;

    end


    % STEP 3: Powells method: make a new direction
    Sc = (x0 - x_orig)/norm(x0 - x_orig);
    for i3 = 1:length(x0) - 1
        sMat(:, i3) = sMat(:, i3 + 1);
    end
    sMat(:, i3 + 1) = Sc;

    % run a case using conjugate direction % bounding phase
    % delta = 0.01;
    s0 = sMat(:, length(x0));

%     [xLEFT, xRIGHT, ~, ~, ~] = boundingPhase(x0, delta, ...
%         objFun, 'direction', s0, 'display', opts);
%     [xLEFT, xRIGHT, ~, ~, ~, funs] = boundingPhase(objFun, x0, delta, ...
%         'direction', s0, 'display', 'none');
    [xLEFT, xRIGHT, ~, ~, ~, funs] = boundingPhase(objFun, x0, delta, ...
        'direction', s0, 'display', 'none', 'LB', LB, 'UB', UB);
    funEvals = funEvals + funs;

    % golden section phase for alpha star
%     [xopt, Fopt, ~] = goldenSection(xLEFT, xRIGHT, objFun);
    [xopt, Fopt, ~, funs] = goldenSection(objFun, xLEFT, xRIGHT);
    funEvals = funEvals + funs;
    aStar = norm(xopt - x0');
    x0 = xopt';

    % STEP 4: Powells method: check alpha criteria and determine going again
    BigAlphaStar = norm(x0 - x_orig);
    dfopt = fopt1 - Fopt;
    if BigAlphaStar < tolerance
        exitFlag = 1;
    elseif iter > iterMax
        exitFlag = 2;
    elseif abs(dfopt) < dFunTol
        exitFlag = 3;
    end
    fopt1 = Fopt;
    
    if dispFlag == 1
%         disp('S0');
%         disp(s0);
%         disp('Xopt');
%         disp(xopt);
%         disp('Fopt');
%         disp(Fopt);
%         disp('aStar');
%         disp(aStar);
        fprintf('%-10d\t %-15d\t %-15d\t %-15d\n', iter, funEvals, Fopt, BigAlphaStar);
    end
%     disp('ITER');
%     disp(iter);
%     disp('funEvals');
%     disp(funEvals);

%     if disps == 1
%         disp('BigAlphaStar');
%         disp(BigAlphaStar);
%         disp('Conj');
%         disp(Sc);
%     end

% %   Dunny why this is here....
%     if iter == 2
%         exitFlag = 1;
%     end
end

end