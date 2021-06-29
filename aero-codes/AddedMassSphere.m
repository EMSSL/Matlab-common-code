function [T, M, k] = AddedMassSphere(a_m, rho, vel)
% Calculate the added mass effect for an oblate spheroid
% Reference: "A review of added mass and fluid inertial forces"
% INPUTS
%     a, b    semimajor axes of the spheroid
%     rho     density, kg/m3
% OUTPUTS
%     T       J, kinetic energy of the fluid motion
%     M       kg, added mass
%     k       
%
% Written by Christopher D. Yoder for the balloon/sail project in 2020. 

% calculate
T = (np.pi/3)*rho*a_m*a_m*a_m*vel*vel;
M = (4/3)*np.pi*rho*a_m*a_m*a_m;
k = 0.5;

end