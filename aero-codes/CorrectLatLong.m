function [dlat, dlong] = CorrectLatLong(lat, long, varargin)
% Function to correct latitude and longitude angles 
% INPUT:
%   lat     = [deg] latitude, scalar or vector
%   long    = [deg] longitude, scalar or vector
% OUTPUT: 
%   dlat    = [deg] latitude, [-90, 90]
%   dlong   = [deg] longitude
% VARARGIN:
%   'long-format'   =   format for longitude degrees
%           'east' (default) bounded between [0, 360], deg East
%           'west' bounded between [-360, 0], deg West
%           'east-west' bounded between [-180, 180], deg East and West
%
% Christopher D. Yoder
% Circa 2018


% error handling
if isnumeric(lat) ~= 1 || isnumeric(long) ~= 1
    error('Inputs are not numeric.');
elseif sum(isnan(lat)) ~= 0 || sum(isnan(long)) ~= 0
    error('Inputs are NaN.');
elseif length(lat) ~= length(long)
    error('Inputs not the same length.');
end

% varargin for formatting
flag = 0;
nanFlag = 0;
if isempty(varargin) == 0
    for i1 = 1:length(varargin)/2
        switch varargin{2*(i1 - 1) + 1}
            case 'long-format'
                
                switch varargin{2*(i1 - 1) + 2}
                    case 'east-west'
                        flag = 1;
                    case 'east'
                        flag = 2;
                    case 'west'
                        flag = 3;
                    otherwise
                        error(['Format ', varargin{2*(i1 - 1) + 2}, ' is not recognized.']);
                end
                
            case 'insert-nan'
                if strcmp(varargin{2*(i1 - 1) + 2}, 'on')
                    nanFlag = 1;
                end
                
            otherwise
                error(['Argument ', varargin{2*(i1 - 1) + 1}, ' is not recognized.']);
        end
    end
end                        
                        

% correct stuff
dlat = NaN*lat;
dlong = NaN*long;
for i1 = 1:length(lat)
    
    % get stuff
    LAT = lat(i1);
    LONG = long(i1);
    
    % correct the latitude 
    if LAT > 90             % if lat is > 90
        LAT = 180 - LAT;    % move the lat to the other side
        LONG = LONG + 180;  % shift the longitude since crossed the world
        
    end
    
    if LAT < -90                    % if lat is less than -90
        LAT = -90 - (LAT + 90);     % shift 
        LONG = LONG + 180;          % correct
    end
    
    % now correct the longitude
    if LONG > 360
        LONG = mod(LONG, 360);
    end
    
    if LONG < 0
        LONG = 360 + LONG;
    end
    
    % assign outputs
    dlat(i1) = LAT;
    dlong(i1) = LONG;
    
end


% post process stuff
switch flag
    case 1     % user wants longitude between [-180, 180] deg
        dlong = dlong - 180;
        
    case 2
        
    case 3
        dlong = dlong - 360;
end


    
end