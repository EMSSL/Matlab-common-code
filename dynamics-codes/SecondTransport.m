function [eqns] = SecondTransport(pos, angvel, t)
% [eqns] = SecondTransport(pos, vel, angvel, angaccl) returns the three-dimensional vector
% corresponding to the inertial derivatives
% pos = B frame position vector
% vel = B frame derivative of pos
% accl = B frame second derivative of pos
% angvel = B frame angular velocity vector
% angaccl = B frame angular acceleration

% eqns = diff(pos, t, t) + 2*cross(angvel, diff(pos, t)) + ...
%     cross(diff(angvel, t), pos) + cross(angvel, cross(angvel, pos));
eqns = simplify(diff(pos, t, t) + 2*CrossMe(angvel, diff(pos, t)) + ...
    CrossMe(diff(angvel, t), pos) + CrossMe(angvel, CrossMe(angvel, pos)));

end