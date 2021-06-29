% Bounding phase algorithm for bounding 1-D optimization
%
% [xLEFT, xRIGHT, fLEFT, fRIGHT, iter, funEvals, mssgFlag] = 
% boundingPhase(objFun, x0, step_size, varargin) returns the left and right 
% bounds (and the required iteration count) for a convex space resulting 
% from the initial guess x0, the step_size, and the desired objective 
% function. fLEET and fRIGHT are the values of the objective function 
% values at xLEFT and xRIGHT.
%   -----------------------------------------------------------------
% INPUTS:
%   objFun      =   objective function for evaluation, must be a string
%   x0          =   starting location
%   delta       =   step size for use in bounding
%   -----------------------------------------------------------------
% OPTIONS:
%   'display'   =   display the information for each iteration in real
%                   time. Default is 'none'. Can be 'iter' for iteration
%                   information.
%   'plot'      =   plot the iterations and display the left and right
%                   values of the function to show bounding. 'on' displays
%                   the plot, 'off' is default.
%   'direction' =   direction of search if this is an n-dimensional
%                   problem. The default is along x1 (1-D search). Specify
%                   as a unit vector in the desired direction.
%   'LB'        =   default is -inf. Lower bounds on the design variables.
%   'UB'        =   default is +inf. Upper bounds on the design variables.
%   -----------------------------------------------------------------
% OUTPUTS:
%   xLEFT       =   left bound 
%   xRIGHT      =   right bound
%   fLEFT       =   function value at left bound
%   fRIGHT      =   function value at right bound
%   iter        =   interations performed
%   funEvals    =   number of function evaluations 
%   mssgFlag    =   exit message flag
%   -----------------------------------------------------------------

function [xLEFT, xRIGHT, fLEFT, fRIGHT, iter, funEvals, mssgFlag] = ...
    boundingPhase(objFun, x0, delta, varargin)
% Written by Christopher D. Yoder 09/15/16
%
% Modifications
% 1. CDY, 7/6/2017, Added capability to handle (poorly probably) the bounds
%    to the design string. Should not exceed the upper and lower bounds when
%    choosing an alphaStar
%
% To Do
% 1. Add exit message flag 

% Initalize and handle inputs
nVars = length(x0);         % get the number of variables for problem
direct = zeros(1, nVars);   % initialize the direction of search
direct(1) = 1;
dispFlag = 0;
plotFlag = 0;
LB = -inf*ones(1, nVars);
UB = inf*ones(1, nVars);
mssgFlag = 'Complete, no errors.';
funEvals = 0;
dFunTol = 1e-6;
if isempty(varargin) ~= 1
    for i = 1:length(varargin)/2
        if strcmp(char(varargin(2*(i - 1) + 1)), 'display') == 1
            if strcmp(char(varargin(2*(i - 1) + 2)), 'iter') == 1
                dispFlag = 1;
            end
        elseif strcmp(char(varargin(2*(i - 1) + 1)), 'plot') == 1
            if strcmp(char(varargin(2*(i - 1) + 2)), 'on') == 1
                plotFlag = 1;
            end
        elseif strcmp(char(varargin(2*(i - 1) + 1)), 'direction') == 1
            direct = cell2mat(varargin(2*(i - 1) + 2));
        elseif strcmp(char(varargin(2*(i - 1) + 1)), 'LB') == 1
            % boundsVector = cell2mat(varargin(2*(i - 1) + 2));
            LB = cell2mat(varargin(2*(i - 1) + 2));
        elseif strcmp(char(varargin(2*(i - 1) + 1)), 'UB') == 1
            % boundsVector = cell2mat(varargin(2*(i - 1) + 2));
            UB = cell2mat(varargin(2*(i - 1) + 2));
            
        else
            error(['Input parameter ', char(varargin(1)), ' is invalid.']);
        end
    end
end
exitFlag = 0;       % initialize exit flag
iter = 1;           % set iter count to one
contFlag = 1;       % initialize the exit flag
modFlag = 0;

if abs(1 - norm(direct)) >= 1e-4
    error(['2 norm of direction vector (', num2str(norm(direct)), ...
        ') ~= 1. Please use unit vectors.']);
end
[m5, n5] = size(direct);
if m5 > n5
    direct = direct';
end
[m5, n5] = size(x0);
if m5 > n5
    x0 = x0';
end

% get first bounds in x 
x0(iter, :) = x0;                               % find the initial x
xL(iter, :) = x0(iter, :) - abs(delta)*direct;  % find the lower x
xR(iter, :) = x0(iter, :) + abs(delta)*direct;  % fine the upper x

% error check for bounding
for i2 = 1:length(x0(1, :))
    if or(x0(iter, i2) < LB(i2),  x0(iter, i2) > UB(i2))
        error('x0 string exceeds the bounds.');
    elseif x0(iter, i2) < LB(i2)
        x0(iter, i2) = LB(i2);
    elseif x0(iter, i2) > UB(i2)
        x0(iter, i2) = UB(i2);
    end
end

% FL(iter) = eval([objFun, '(xL(iter, :))']);  % find the initial f
% F0(iter) = eval([objFun, '(x0(iter, :))']);  % find the lower f
% FR(iter) = eval([objFun, '(xR(iter, :))']);  % find the upper f


FL(iter) = objFun(xL(iter, :));
F0(iter) = objFun(x0(iter, :));
FR(iter) = objFun(xR(iter, :));
funEvals = funEvals + 3;

if FL(iter) >= F0(iter) && F0(iter) >= FR(iter) % if descenting left to right
    delta = abs(delta);                         % delta is positive
elseif FR(iter) >= F0(iter) && F0(iter) >= FL(iter) % if ascending left to right
    delta = -delta;                             % delta is negative
else
    xLEFT = xL(iter, :);         % otherwise you win
    xRIGHT = xR(iter, :);        % set up the bounds
    fLEFT = FL(iter);            % and report
    fRIGHT = FR(iter);           % the values if asked
    contFlag = 0;                % and exit
    iter = iter - 1;
end

if contFlag == 1                                                % if you didn't win
    x0(iter + 1, :) = x0(iter, :) + (2^iter)*delta*direct;      % calculate new x value
    % F0(iter + 1) = eval([objFun, '(x0(iter + 1, :))']);       	% find new f value
    
    % error check for bounding
    for i2 = 1:length(x0(1, :))
        if x0(iter + 1, i2) < LB(i2)
            x0(iter + 1, i2) = LB(i2);
        elseif x0(iter + 1, i2) > UB(i2)
            x0(iter + 1, i2) = UB(i2);
        end
    end
    
    F0(iter + 1) = objFun(x0(iter + 1, :));
    funEvals = funEvals + 1;
    
    while exitFlag ~= 1
            
        if F0(iter + 1) + dFunTol < F0(iter)  % if the new values are less than current value
            iter = iter + 1;        % you're still going down hill buddy
            x0(iter + 1, :) = x0(iter, :) + (2^iter)*delta*direct; % find a new x
            % F0(iter + 1) = eval([objFun, '(x0(iter + 1, :))']);    % find a new f
            
            % error check for bounding
            for i2 = 1:length(x0(1, :))
                if x0(iter + 1, i2) < LB(i2)
                    x0(iter + 1, i2) = LB(i2);
                elseif x0(iter + 1, i2) > UB(i2)
                    x0(iter + 1, i2) = UB(i2);
                end
            end
            
            F0(iter + 1) = objFun(x0(iter + 1, :));
            funEvals = funEvals + 1;
        else
            exitFlag = 1;
        end      
    end

    % assign values
    if iter == 1
        xLEFT = x0(iter, :);
        xRIGHT = x0(iter + 1, :);
        fLEFT = F0(iter);
        fRIGHT = F0(iter + 1);
        % disp('Modification enabled');
        mssgFlag = 'Modification enabled.';
    else
        xLEFT = x0(iter - 1, :);       % label things for the outputs
        xRIGHT = x0(iter + 1, :);
        fLEFT = F0(iter - 1);
        fRIGHT = F0(iter + 1);
    end
end


% plot to check things
if plotFlag == 1
    plot(x0, F0, '-*')
    grid on
    hold on
    plot(xLEFT, fLEFT, 'o', xRIGHT, fRIGHT, 'o')
end

% display iteration table
if dispFlag == 1
    ITER = 1:1:iter + 1;
    ITER = ITER';
    % x0 = x0';
    F0 = F0';
    T = table(ITER, x0, F0);
    disp('Bounding Phase Results');
    disp(T);
end


end