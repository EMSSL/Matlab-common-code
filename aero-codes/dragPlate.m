% Fucntion to return the drag coefficient for flow normal to a flat plate. 
% Written by Christopher D. Yoder for EMSSL 05/21/2016.
%
% INPUTS:
%   Rynlds      =   Reynolds number of the flow normal to the plate
%   
% OUTPUTS:
%   Cd          =   drag coefficient of the plate
%
% REFERENCE:
%   [1] Hoerner, S. F., Fluid-Dynamic Drag.pdf, 1965.
% 

function [Cd] = dragPlate(Rynlds)

% error handling
if isnumeric(Rynlds) ~= 1
    error(['Reynolds number is not numeric. Please pass in a ', ...
        'number to this function.']);
end

% tablS = [Re, Cd], derived from Hoerner drag book 
tablS = [0.01       2300
         0.1        225
         1          23
         3          8.5
         5          5.5
         8          4.3
         10         3.6
         20         2.45
         50         1.7
         80         1.5
         100        1.5
         200        1.75
         300        2
         450        1.5
         800        1.25
         1000       1.45
         2000       1.2
         3500       1.2
         6000       1.15
         10000      1.17
         100000     1.17
         1000000	1.17];   
   
% interpolate using the table above to find the value of Cd     
Cd = interp1(tablS(:, 1), tablS(:, 2), Rynlds, 'pchip');
end