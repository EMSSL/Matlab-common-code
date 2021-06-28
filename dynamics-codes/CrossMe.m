function [prod] = CrossMe(vec1, vec2)
% [prod] = CrossMe(vec1, vec2) returns the cross product for two vectors
%        | i   j   k|
% vec1 = [x1, y1, z1]
% vec2 = [x2, y2, z2]
% thus:
%  vec1 x vec2 = [ y1*z2 - y2*z1] = [ vec1(2)*vec2(3) - vec2(2)*vec1(3)]
%                [-x1*z2 + x2*z1] = [-vec1(1)*vec2(3) + vec2(1)*vec1(3)]
%                [ x1*y2 - x2*y1] = [ vec1(1)*vec2(2) - vec2(1)*vec1(2)]
%
% This is because the cross() function in MATLAB no worky for symbolics
% 
% INPUTS
%   vec1, vec2      = inputs for vec1 x vec2
% OUTPUTS
%   prod            = resulting cross product vector

prod = sym('prod', [3, 1]);
prod(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2);
prod(2) = -vec1(1)*vec2(3) + vec1(3)*vec2(1);
prod(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1);

end
