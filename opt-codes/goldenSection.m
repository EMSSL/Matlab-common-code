% Function to perform 1-D optimization using Golden section method.
%  
% goldenSection(objFun,xL,xR) returns the value of the optimum point of the function
% objFun when bounded between left bound xL and right bound xR. 
% 
% [xopt, Fopt, iter, funEvals] = goldenSection(objFunxL, xL, xR) returns the value of the optimum point xopt 
% of the function objFun when bounded between left bound xL and right bound xR
% and the corresponding value of the objective function. iter lists the
% number of iterations used to get to the solution.
%
% [xopt, Fopt, iter, funEvals] = goldenSection(objFun,xL,xR,...) also runs the golden
% section method subjected to the options specified. Options are given
% below. 
%   -----------------------------------------------------------------
% INPUTS:
%     objFun  =   Objective function to evaluate. Must be a function handle
%     xL      =   Left bound from bounding phase
%     xR      =   Right bound from bounding phase
%   -----------------------------------------------------------------
% OPTIONS:
%     'display'   
%             =   display the information for each iteration in real
%                 time. Default is 'none'. Can be 'iter' for iteration
%                 information.
%     'term'  =   Termination criteria for golden section. If not specified,
%                 the default is 1e-5
%     'direction' 
%             =   direction of search if this is an n-dimensional
%                 problem. The default is along x1 (1-D search). Specify
%                 as a unit vector in the desired direction.
%   -----------------------------------------------------------------
% OUTPUTS:
%     xopt      =   output solution for the problem
%     Fopt      =   value of the objective function at the point xopt
%     iter      =   number of iterations (function calls) to reach the
%                   solution
%     funEvals  =   number of function evaluations
%   -----------------------------------------------------------------

function [xopt, Fopt, iter, funEvals] = goldenSection(objFun, xL, xR, varargin)
% Written by S. Ferguson for MAE 531. 
% Modified 9/15/16 by Christopher D. Yoder
% 
% Modifications:
% 1. CDY, 9/15/17, Added varargin syntax for display usage
% 2. CDY, 6/30/17, Modified input order to match other functions, modified
%    the objFun to accept function handles, not strings.
% 3. CDY, 6/30/17, Added function evaluations to outputs.
% 
% To Do:
% 1. Change the display to display as the function runs, not afterwards.

% Set up inputs
ex_terminate = 0.00001;     % default termination criteria
n5 = length(xL);
direct = zeros(1, n5);
direct(1) = 1;
dispFlag = 0;
funEvals = 0;
if isempty(varargin) ~= 1            % if user chooses termination criteria
    for i = 1:length(varargin)/2
        if strcmp(char(varargin(2*(i - 1) + 1)), 'term') == 1   % and if specified
            ex_terminate = cell2mat(varargin(2*(i - 1) + 2)); % save criteris
        elseif strcmp(char(varargin(2*(i - 1) + 1)), 'direction') == 1
            direct = cell2mat(varargin(2*(i - 1) + 2));
        elseif strcmp(char(varargin(2*(i - 1) + 1)), 'display') == 1
            if strcmp(char(varargin(2*(i - 1) + 2)), 'iter') == 1
                dispFlag = 1;
            end
        else
            % else throw an error
            error(['Input parameter ', char(varargin(2*(i - 1) + 1)),...
                ' is invalid.']);
        end
    end
end

% Generic golden section code
iter = 1;
xR(iter, :) = xR;
xL(iter, :) = xL;
% FR(iter) = eval([objFun, '(xR(iter, :))']);
% FL(iter) = eval([objFun, '(xL(iter, :))']);
FR(iter) = objFun(xR(iter, :));
FL(iter) = objFun(xL(iter, :));
L_original = norm(xR(iter, :) - xL(iter, :));       % initial L value, save for later
L(iter) = L_original;             % initialize L for while loop
tau = 0.381966;             % constant, doesn't change
funEvals = funEvals + 2;

x1(iter, :) = (1-tau)*xL(iter, :) + tau*xR(iter, :);       % calculate x1
% F1(iter) = eval([objFun, '(x1(iter, :))']);    % calculate function at x1
F1(iter) = objFun(x1(iter, :));    % calculate function at x1
x2(iter, :) = tau*xL(iter, :) + (1-tau)*xR(iter, :);       % calulate x2
% F2(iter) = eval([objFun, '(x2(iter, :))']);    % calculate function at x2
F2(iter) = objFun(x2(iter, :));    % calculate function at x1
funEvals = funEvals + 2;

iter = iter + 1;                   % set iter to two
convergence_ratio = 1;      % set convergence ratio

% Ready, ...set, ...optimise!
while (convergence_ratio > ex_terminate)    
    
    if (F1(iter - 1) > F2(iter - 1))    % if f1 is greater
        xL(iter, :) = x1(iter - 1, :);    % eliminate f1 side
        xR(iter, :) = xR(iter - 1, :);
        FL(iter) = F1(iter - 1);    % rename variables
        FR(iter) = FR(iter - 1);
        x1(iter, :) = x2(iter - 1, :);
        F1(iter) = F2(iter - 1);
        x2(iter, :) = tau*xL(iter, :) + (1-tau)*xR(iter, :);
        % F2(iter) = eval([objFun, '(x2(iter, :))']);
        F2(iter) = objFun(x2(iter, :));    % calculate function at x1
        iter = iter + 1;
        funEvals = funEvals + 1;
       
    else
        xR(iter, :) = x2(iter - 1, :);    % if f2 is greater
        xL(iter, :) = xL(iter - 1, :);
        FR(iter) = F2(iter - 1);
        FL(iter) = FL(iter - 1);
        x2(iter, :) = x1(iter - 1, :);
        F2(iter) = F1(iter - 1);
        x1(iter, :) = (1-tau)*xL(iter, :) + tau*xR(iter, :);
        % F1(iter) = eval([objFun, '(x1(iter, :))']);
        F1(iter) = objFun(x1(iter, :));    % calculate function at x1
        iter = iter + 1;
        funEvals = funEvals + 1;
  
    end   

%     xL
%     FL
%     x1
%     F1
%     x2
%     F2
%     xR
%     FR
    
    L(iter - 1) = norm(xR(iter - 1, :) - xL(iter - 1, :));
    convergence_ratio = L/L_original;
end

% outputs
xopt = (xL(iter - 1, :) + xR(iter - 1, :))/2;
% Fopt = eval([objFun, '(xopt)']);
Fopt = objFun(xopt);
funEvals = funEvals + 1;

% display table stuff
if dispFlag == 1
    ITER(:, 1) = 1:1:iter - 1;
    F1 = F1';
    F2 = F2';
    T = table(ITER, x1, x2, F1, F2);
    disp('Golden Section Results');
    disp(T);
end

end