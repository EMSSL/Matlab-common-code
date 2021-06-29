% Fucntion to return the drag coefficient for flow normal to a sphere. 
% Written by Christopher D. Yoder for EMSSL 05/21/2016.
%
% INPUTS:
%   Rynlds      =   Reynolds number of the flow normal to the sphere
%   
% OUTPUTS:
%   Cd          =   drag coefficient of the sphere
%
% REFERENCE:
%   [1] Hoerner, S. F., Fluid-Dynamic Drag.pdf, 1965.


function [Cd] = dragSphe(Rynlds)

% error handling
if isnumeric(Rynlds) ~= 1
    error(['Reynolds number is not numeric. Please pass in a ', ...
        'number to this function.']);
end

% tablC = [Re, Cd]
tablS = [0.019      1000
         0.03       650
         0.05       400
         0.1        190
         0.5        45
         1          25
         2          14
         5          7.25
         10         4.5
         20         2.8
         40         1.8
         60         1.45
         80         1.3
         250        0.75
         1000       0.475
         2200       0.41
         20000      0.5
         50000      0.5
         100000     0.5
         200000     0.5
         325000     0.5
         500000     0.1
         1000000	0.1
         1e9        0.1];   
     
% Cd = interp1(tablS(:, 1), tablS(:, 2), Rynlds, 'pchip');
% Cd = interp1(tablS(:, 1), tablS(:, 2), Rynlds, 'spline');
Cd = interp1(tablS(:, 1), tablS(:, 2), Rynlds, 'nearest', 'extrap');

end