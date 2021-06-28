function [vec] = MakeVec(x, y, z)
% [vec] = MakeVec(x, y, z) returns a vector vec of elements x, y, z. 

% error handle, check type
if isa([x, y, z], 'sym') == 1
    vec = sym('vecz', [3, 1]);
elseif isnumeric(x) == 1
    vec = NaN(3, 1);
else
    error('x neither symbol nor numeric.');
end

% populate
vec(1) = x;
vec(2) = y;
vec(3) = z;

end