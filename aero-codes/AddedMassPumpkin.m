function [T_v, M_v, k_v, T_h, M_h, k_h] = AddedMassPumpkin(a_m, b_m, rho, vel_v, vel_h)
% Calculate the added mass effect for an oblate spheroid
% Reference: Brennen, CE. A review of added mass and fluid inertial forces
%           1982, https://apps.dtic.mil/sti/citations/ADA110190
% INPUTS
%     a, b      semimajor axes of the spheroid, a>b (e.g. b=hgt/2, a=diam/2)
%     rho       density, kg/m3
%     vel_v, _h velocity, m/s, in vertical and horizontal
% OUTPUTS
%     T         J, kinetic energy of the fluid motion
%     M         kg, added mass
%     k         ratio?
% 
% Written by Christopher D. Yoder for the balloon/sail project in 2020. 


% check a < b 
if a_m < b_m
    error("a should be greater than b.")
end

% correct notation from previous literature
b_m1 = a_m;
a_m1 = b_m;
a_m = a_m1;
b_m = b_m1;

% vertical motion, k1
e = sqrt(1 - b_m*b_m/(a_m*a_m));
if e == 0
    e = 1e-12;
end
alpha0 = (e - sqrt(1 - e*e)*asin(e))*(2/(e*e*e));
beta0 = (sqrt(1 - e*e)*asin(e) - e*(1 - e*e))*(1/(e*e*e));
T_v = (2/3)*rho*pi*a_m*b_m*b_m*vel_v*vel_v*alpha0/(2 - alpha0);
M_v = (4/3)*rho*pi*a_m*b_m*b_m;
k_v = alpha0/(2 - alpha0);
% print(k_v)

% horizontal motion, k2
e = sqrt(1 - b_m*b_m/(a_m*a_m));
if e == 0
    e = 1e-12;
end
alpha0 = (e - sqrt(1 - e*e)*asin(e))*(2/(e*e*e));
beta0 = (sqrt(1 - e*e)*asin(e) - e*(1 - e*e))*(1/(e*e*e));
T_h = (2/3)*rho*pi*a_m*b_m*b_m*vel_h*vel_h*beta0/(2 - beta0);
M_h = (4/3)*rho*pi*a_m*b_m*b_m;
k_h = beta0/(2 - beta0);

end