% [Vx, Vy, Vz] = downWash(point, vStart, vEnd) computes the
% values of velocity (downWash) {Vx; Vy; Vz} for a point with coordinates 
% {x; y; z} subjected to a vortex segment with starting point (vStart) 
% {x1; y1; z1} and end point (vEnd) {x2; y2; z2}. The starting point and
% ending point are used to denote the direction of circulation of the vortex.
% 
% TYPE 1: finite vortex, inputs are point, vStart, vEnd coordinates
% TYPE 2: semiinfinite vortex with gamma along t, inputs are point, vStart, trail direction vector  
% TYPE 3: semiinfinite vortex with gamma against t, inputs are point, vStart, trail direction vector
% For semi-infinite vortex filaments, the vEnd vector is replaced with Inf.
% An additional direction argument is the direction vector of 
% circulation. THIS DIRECTION DENOTES THE CURL!!!

% Written by Christopher D. Yoder for MAE 561: Wing Theory Homework
% assignment one, problem two.
%
% References: 
% Course material provided for vector form of the Biot Savart Law.



function [Vx, Vy, Vz] = downWash1(type, varargin) 

if type == 1
    % TYPE 1: varargin = [point], [vStart], [vEnd]
    point = cell2mat(varargin(1));
    vStart = cell2mat(varargin(2));
    vEnd = cell2mat(varargin(3));
    a = vStart - point;
    b = vEnd - point;
    axb = cross(a, b);
    product = axb*(norm(a) + norm(b))*(1 - ...
        dot(a, b)/(norm(a)*norm(b)))/dot(axb, axb);
    
    % catch points which are colinear
    if norm(axb) <= 1e-4
        product = zeros(1, 3);
    end

elseif type == 2
    % TYPE 2: varargin = [point], [vStart], [trailDirectionVector]
    point = cell2mat(varargin(1));
    vStart = cell2mat(varargin(2));
    t = cell2mat(varargin(3))/norm(cell2mat(varargin(3)));
    a = vStart - point;
    axt = cross(a, t);
    product = axt*(1 - dot(a, t)/norm(a))/(dot(axt, axt));

    % catch points which are colinear
    if norm(axt) <= 1e-4
        product = zeros(1, 3);
    end
    
elseif type == 3
    % TYPE 3: varargin = [point], [vStart], [trailDirectionVector]
    point = cell2mat(varargin(1));
    vStart = cell2mat(varargin(2));
    t = cell2mat(varargin(3))/norm(cell2mat(varargin(3)));
    a = vStart - point;
    axt = cross(a, t);
    product = -(axt*(1 - dot(a, t)/norm(a))/(dot(axt, axt)));
    
    % catch points which are colinear
    if norm(axt) <= 1e-4
        product = zeros(1, 3);
    end    
    
else
    error(['Type ', num2str(type), ' is not defined.']);
end

% assign outputs
Vx = product(1)/(4*pi);
Vy = product(2)/(4*pi);
Vz = product(3)/(4*pi);

end