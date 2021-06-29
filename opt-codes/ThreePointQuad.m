% Three-point-quadrature code for 1-D optimization
%
% [xOPT, fOPT, convergence, iter] = ThreePointQuad(xLEFT, xRIGHT, objFun, varargin)
% returns the optimal value of x, xOPT, the optimized function value, fOPT,
% the convergence criteria as satisfied, and the number of iterations used.
% xLEFT and xRIGHT are the x values which bound the space taken from boundingPhase.m
% of the objective function objFun. 
%
% INPUTS:
%   xLEFT       =   left location in x
%   xRIGHT      =   right location in x
%   objFun      =   objective function for evaluation, must be a string
%
% OPTIONS:
%   'tolerance' =   convergence criteria tolerance. Default is 1e-4. 
%   'display'   =   return information used during the optimization.
%                   Default is 'none'. 'iter' displays the information 
%                   for each iteration. 
%   'plotIter'  =   Display a plot of the three points and the quadratic
%                   fit for a given iteration increment value. Default is
%                   'off'. Followed by 'i' where i indicates plots be made on
%                   every ith iteration. 
%
% OUTPUTS:
%   xOPT        =   optimum x value 
%   fOPT        =   optimum function value
%   convergence =   convergence achieved
%   iter        =   interations performed

function [xOPT, fOPT, convergence, iter] = ThreePointQuad(xLEFT, xRIGHT, objFun, varargin)
% Written by Christopher D. Yoder for MAE 531 EDO

% Initalize and handle inputs
dispFlag = 0;
epx = 1e-4;
plotFlag = 0;
xL0 = xLEFT;
xR0 = xRIGHT;
if isempty(varargin) ~= 1
    for i = 1:length(varargin)/2
        if strcmp(char(varargin(2*(i - 1) + 1)), 'display') == 1
            if strcmp(char(varargin(2*(i - 1) + 2)), 'iter') == 1
                dispFlag = 1;
            end
        elseif strcmp(char(varargin(2*(i - 1) + 1)), 'tolerance') == 1
            epx = cell2mat(varargin(2*(i - 1) + 2));
        elseif strcmp(char(varargin(2*(i - 1) + 1)), 'plotIter') == 1
            plotFlag = cell2mat(varargin(2*(i - 1) + 2));
        else
            error(['Input parameter ', char(varargin(1)), ' is invalid.']);
        end
    end
end

% Get defaults
xMID = (xLEFT + xRIGHT)/2;
fLEFT = eval([objFun, '(xLEFT)']);
fMID = eval([objFun, '(xMID)']);
fRIGHT = eval([objFun, '(xRIGHT)']);
exitFlag = 0;
iter = 1;

% ready...set...optimize!
while exitFlag ~= 1
    
% Find minimums
xVEC(iter, :) = [xLEFT, xMID, xRIGHT];
fVEC(iter, :) = [fLEFT, fMID, fRIGHT];
[fMIN, minInd] = min(fVEC(iter, :));
xMIN = xVEC(iter, minInd);

% Find constants
a0(iter, 1) = fLEFT;
a1(iter, 1) = (fMID - fLEFT)/(xMID - xLEFT);
a2(iter, 1) = (1/(xRIGHT - xMID))*((fRIGHT - fLEFT)/(xMID - xLEFT) - ...
    (fMID - fLEFT)/(xMID - xLEFT));
xSTAR(iter, 1) = (xMID + xLEFT)/2 - a1(iter, 1)/(2*a2(iter, 1));
fSTAR(iter, 1) = eval([objFun, '(xSTAR(iter, 1))']);

% check convergence criteria
if abs((xSTAR(iter, 1) - xMIN)/xMIN) <= epx
    exitFlag = 1;
end

% Save best values
if xMIN < xSTAR(iter, 1)
    xVEC1 = [xLEFT, xMIN, xSTAR(iter, 1), xRIGHT];
    fVEC1 = [fLEFT, fMIN, fSTAR(iter, 1), fRIGHT];
else
    xVEC1 = [xLEFT, xSTAR(iter, 1), xMIN, xRIGHT];
    fVEC1 = [fLEFT, fSTAR(iter, 1), fMIN, fRIGHT];
end
[fMID, minInd] = min(fVEC1);
xMID = xVEC1(minInd);
xLEFT = xVEC1(minInd - 1);
xRIGHT = xVEC1(minInd + 1);
fLEFT = fVEC1(minInd - 1);
fRIGHT = fVEC1(minInd + 1);

% Infinite loop protection
iter = iter + 1;
if iter >= 1e5
    break
end

end

% assign outputs
xOPT = xMID;
fOPT = fMID;
convergence = abs((xSTAR(iter - 1, 1) - xMIN)/xMIN);
iter = iter - 1;

% display iteration information
if dispFlag == 1
    ITER(:, 1) = 1:1:iter;
    x1(:, 1) = xVEC(:, 1);
    x2(:, 1) = xVEC(:, 2);
    x3(:, 1) = xVEC(:, 3);
    F1(:, 1) = fVEC(:, 1);
    F2(:, 1) = fVEC(:, 2);
    F3(:, 1) = fVEC(:, 3);
    table(ITER, a0, a1, a2, x1, F1, x2, F2, x3, F3, xSTAR, fSTAR)
end

% plot the iterations requested
if plotFlag ~= 0
    xv = linspace(xL0, xR0, 100);
    % xv = linspace(-5, 5, 100);
    for j = 1:plotFlag:iter + 1
        if j ~= 1
            k = j - 1;
        else 
            k = j;
        end
        figure
        F = eval([objFun, '(xv)']);
        plot(xv, F, '-')
        grid on
        xlabel('x')
        ylabel('F(x)')
        lims = axis;
        hold on
        Fn = a0(k, 1) + a1(k, 1).*(xv - xVEC(k, 1)) + ...
            a2(k, 1).*(xv - xVEC(k, 1)).*(xv - xVEC(k, 2));
        plot(xv, Fn, '-.')
        legend('Function', 'Approximation')
        title(['Iteration ', num2str(k)])
        axis(lims)
        hold off
    end
end
end