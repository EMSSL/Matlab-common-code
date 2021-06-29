function testCode()
% test code for pattern search algorithm

clc
close all
clear


% 2D example with plotting
% grid solution
x1 = -5:0.25:5;
y1 = -5:0.25:5;
for i1 = 1:length(x1)
    for j1 = 1:length(y1)
        z1(i1, j1) = obj2D([x1(i1), y1(j1)]);
    end
end
contour(x1, y1, z1, linspace(-12, 120, 39))
% contourcbar;
% cbar
% set(gcf, 'Position', [1241         462         649         493]);

% PowellsMethod search
x0 = [1, -1];
[xopt, Fopt, iter, funEvals] = PowellsMethod(@(x)obj2D(x), x0);
xopt
hold on
plot(xopt(1), xopt(2), 'rx', 'MarkerSize', 10)

[xopt, Fopt, iter, funEvals, exitFlag] = PatternSearch(@(x)obj2D(x), x0)
xopt
plot(xopt(1), xopt(2), 'bo', 'MarkerSize', 10)

end

function F = obj2D(x)
% function for evaluation

% inputs
x1 = x(1);
x2 = x(2);
x3 = -1;

% OBJECTIVE FUNCTION
F = x1*x1 + 2*x2*x2 + 2*x3*x3 + 2*x1*x2 + 2*x2*x3;
end