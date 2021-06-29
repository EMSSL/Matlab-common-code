function [xopt, Fopt, iter, funEvals, exitFlag] = PatternSearch(objFun, x0, varargin)
% [xopt, Fopt, iter, funEvals] = PatternSearch(objFun, x0, varargin)
% returns the optimal design string, the function evaluation, number of
% iterations, and function evaluations for the objecive function and design
% string using the pattern search algorithm. 
%
%

% Written by Christopher D. Yoder on 7/6/2017

% Handle inputs and varargin
delta = 0.1*ones(1, length(x0));    % step size in each direction
if isrow(x0) ~= 1
    x0 = xo';
end
m1 = length(x0);
contFlag = 1;                   % 1 = continuous domain, 0 = discrete domain
exitFlag = 0;                   % 1 = kickout of the loop
tol = 1e-3;                     % function tolerance
iterMax = 1e3;
minDelta = 1e-2;                % minimum delta allowed

% STEP 1: Evaluate starting point 
F = [];
x = [];
iter = 1;
x(iter, :) = x0;                % get first design string
funEvals = 0;
F(iter) = objFun(x(iter, :));   % get function value
funEvals = funEvals + 1;
minVec = NaN(1, length(x0));
s1 = eye(length(x0));

% find other positions
while exitFlag == 0
    
    x(iter, :)
      
    % get surrounding points
    xP = x(iter, :) + delta;
    xM = x(iter, :) - delta;
    FP = NaN(1, length(x0));
    FM = NaN(1, length(x0));
    deltaF = zeros(1, length(x0));
    
    % plot - testCase
    hold on
    plot(x(iter, 1), x(iter, 2), 'ko')
    plot([xM(1), xP(1)], x(iter, 2)*[1, 1], 'k-o')
    plot(x(iter, 1)*[1, 1], [xM(2), xP(2)], 'k-o')
    pause(1)
    minG = 1;
    cntr = 0;
    
    for i1 = 1:m1                   % for all design variables
        
        % get positive increment
        FP(i1) = objFun(x(iter, :) + delta.*s1(i1, :));
        funEvals = funEvals + 1;
        
        % get negative increment
        FM(i1) = objFun(x(iter, :) - delta.*s1(i1, :));
        funEvals = funEvals + 1;
        
        % choose which step to go with
        if and(F(iter) < FP(i1) - tol, F(iter) < FM(i1) - tol) ~= 1      % when the current is greater than others....
            if FP(i1) <= FM(i1)                 % if FP gives better function value
                minL = abs(FP(i1)) - abs(FM(i1));
                if minL < minG
                    minG = minL;
                    cntr = i1;
                end
            elseif FP(i1) > FM(i1)              % if FM gives better function value
                minL = abs(FM(i1)) - abs(FP(i1));
                if minL < minG
                    minG = minL;
                    cntr = -i1;
                end
            end     
        end
    end
    
    
    % check exit conditions - update and evaluate
    if cntr == 0
        delta = delta/2;
        x(iter + 1, :) = x(iter, :);
    else
        x(iter + 1, :) = x(iter, :);
        if cntr > 0
            x(iter + 1, cntr) = xP(cntr);
        else
            x(iter + 1, -cntr) = xM(-cntr);
        end
    end
    F(iter + 1) = objFun(x(iter + 1, :));
    funEvals = funEvals + 1;
    
    % exit flag
    if F(iter + 1) > F(iter)
        exitFlag = 1;
    elseif iter > iterMax
        exitFlag = 2;  
    elseif delta < minDelta
        exitFlag = 3;
    end
    iter = iter + 1;
    
    
end

% assign outputs
xopt = x(end, :);
Fopt = F(end);

end