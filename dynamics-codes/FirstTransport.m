function [eqns] = FirstTransport(pos, ang, t)
% [eqns] = FirstTransport(pos, ang) returns the three-dimensional vector
% corresponding to the inertial derivatives
% pos = B frame position vector
% vel = B frame derivative of pos
% ang = B frame angular velocity vector

% eqns = simplify(diff(pos, t) + cross(ang, pos));
eqns = simplify(diff(pos, t) + CrossMe(ang, pos));

end