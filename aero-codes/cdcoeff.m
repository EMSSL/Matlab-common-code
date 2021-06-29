function [Cl,Cd]=cdcoeff(rho,V,D,mu,alfa)
% CDCOEFF produces the coefficient of drag for cylinders
% RHO is density in kg/m^3
% V is relative fluid velocity in m/s
% D is reference length in m
% mu is dynamic viscosity in N.s/m^2
% ALFA is angle of attack in degrees
% CL is coefficient of lift 
% CD is coefficient of drag

% angle alpha is defined as:
%      \\ 
%  -->  \\
%  -->   \\
%  FLOW   \\  PLATE
%  -->    /\\
%  -->   /  \\
%      ALPHA \\
%  ____/_____ \\     
% alpha is the small angle between flow and plate inclined

% Calculate Reynolds Number
Re=(rho*V*D)/mu;

% tablC = [Re, Cdvert]
tablC = [0.01       380
         0.025      175
         0.05       98
         0.1        58
         0.5        19
         1          12
         2          7.5
         5          4
         10         2.8
         20         2.1
         40         1.8
         60         1.6
         80         1.5
         100        1.4
         1000       1.08
         4000       1
         10000      1.15
         50000      1.22
         100000     1.23
         200000     1.2
         300000     1
         500000     0.25
         1000000	0.325
         1e16       0.325];    
% % interpolate and calculate the coefficients      
if Re < tablC(1, 1)
    error('Re < min(CD)');
end
if Re > tablC(end, 1)
    error('Re > max(CD)');
end
% find nominal value of vertical tether drag
CdN = interp1(tablC(:, 1), tablC(:, 2), Re);
% correct for incline and act in flow direction
Cd = abs(CdN*((cosd(alfa))^3)); 
% correct for incline and act y = +ve INF
Cl = CdN*((sind(alfa)))*(cosd(alfa)^2);                       
end