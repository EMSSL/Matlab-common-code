function [vec] = MakeMat(xx, yy, zz, xy, yz, xz)
% [vec] = MakeMat(xx, yy, zz, xy, yz, xz) returns a matrix
% vec of elements similar to an inertia tensor. 

% error handle, check type
if isa(xx, 'sym') == 1
    vec = sym('vecz', [3, 3]);
elseif isnumeric(xx) == 1
    vec = NaN(3, 3);
else
    error('x neither symbol nor numeric.');
end

% populate
% IBB = [Ixx, Ixy, Ixz]
%       [Ixy, Iyy, Iyz]
%       [Ixz, Iyz, Izz]
vec(1, 1) = xx;
vec(1, 2) = xy;
vec(1, 3) = xz;
vec(2, 1) = xy;
vec(2, 2) = yy;
vec(2, 3) = yz;
vec(3, 1) = xz;
vec(3, 2) = yz;
vec(3, 3) = zz;

end