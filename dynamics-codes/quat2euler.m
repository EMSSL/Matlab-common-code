function [theta_1, theta_2, theta_3] = quat2euler(q0, q1, q2, q3, choice)
% function to give out euler angle set for a given set of quaternion Parameters  
% INPUTS 
%   q0, q1, q2, q3      = quaternions (twist, i, j, k)
%   choice              = desired mazzoleni rotation sequence
% OUTPUTS
%   theta_1, _2, _3     = euler angles, in rads
% 
% NOTE: 
% Please check if the quaternion parameters are valid !!! (Magnitude < =1 )
%   By SDIVI
% 
% 1 - XYZ   5 - YXZ   9 - ZXY 
% 2 - XZY   6 - YZX  10 - ZYX
% 3 - XYX   7 - YXY  11 - ZXZ
% 4 - XZX   8 - YZY  12 - ZYZ

% Modification table
% CDY, 09-24-2019, Fixed BCO matrix computation

% check symbolic or numeric
keyboard
if strcmp(class(q0),'sym') == 1         % check for symbolic variables
    BCO = sym('BCO',[3,3]);             
elseif isnumeric(q0) == 1               % check for numeric variables
    BCO = zeros(3, 3);
    
    % normalize the quaterions (just in case)
    if sqrt(q0^2 + q1^2 + q2^2 + q3^2) ~= 1
        [q0t,q1t,q2t,q3t] = quat_normalize(q0,q1,q2,q3);
    else
        q0t = q0;
        q1t = q1;
        q2t = q2;
        q3t = q3;
    end
    
else
    warning('q0 is neither symbolic or numeric. Attempting to execute anyway.');
end

% % form the BCO matrix
% % BCO(1,1) = q0t^2 + q1t^2 + q2t^2 + q3t^2 ;
% BCO(1,1) = q0t^2 + q1t^2 - q2t^2 + q3t^2 ;
% BCO(1,2) = 2*(q1t*q2t + q0t*q3t) ;
% BCO(1,3) = 2*(q1t*q3t - q0t*q2t) ;
% BCO(2,1) = 2*(q1t*q2t - q0t*q3t) ;
% % BCO(2,2) = (q0t^2 - q1t^2) + (q2t^2 + q3t^2) ;
% BCO(2,2) = q0t^2 - q1t^2 + q2t^2 - q3t^2 ;
% BCO(2,3) = 2*(q2t*q3t + q0t*q1t);
% BCO(3,1) = 2*(q1t*q3t + q0t*q2t);
% BCO(3,2) = 2*(q2t*q3t - q0t*q1t);
% BCO(3,3) = q0t^2 - q1t^2 - q2t^2 + q3t^2 ;

% % form the BCO matrix - Mazz notation
% BCO(1,1) = q0t^2 + q1t^2 - q2t^2 - q3t^2 ;
% BCO(1,2) = 2*(q1t*q2t + q0t*q3t) ;
% BCO(1,3) = 2*(q1t*q3t - q0t*q2t) ;
% BCO(2,1) = 2*(q1t*q2t - q0t*q3t) ;
% BCO(2,2) = q0t^2 - q1t^2 + q2t^2 - q3t^2 ;
% BCO(2,3) = 2*(q2t*q3t + q0t*q1t);
% BCO(3,1) = 2*(q1t*q3t + q0t*q2t);
% BCO(3,2) = 2*(q2t*q3t - q0t*q1t);
% BCO(3,3) = q0t^2 - q1t^2 - q2t^2 + q3t^2 ;
% % form the OCB matrix
% M = transpose(BCO);

% form BCO matrix - Mazz notation, always use one function!
% M = transpose(BCO) = OCB
BCO = MakeBCO_quat(q0t, q1t, q2t, q3t);
OCB = transpose(BCO);



% compute the Euler angles based on the desired sequence
% WHERE DOES THIS COME FROM????
switch choice
    case 1 
        theta_1 = atan2(-OCB(2,3),OCB(3,3)) ; 
        theta_2 = atan2(OCB(1,3),sqrt(1-OCB(1,3)^2)) ;
        theta_3 = atan2(-OCB(1,2),OCB(1,1)) ;
    case 2 
        theta_1 = atan2(OCB(3,2),OCB(2,2)) ; 
        theta_2 = atan2(-OCB(1,2),sqrt(1-OCB(1,2)^2)) ;
        theta_3 = atan2(OCB(1,3),OCB(1,1)) ;
    case 3 
        theta_1 = atan2(OCB(2,1),-OCB(3,1)) ; 
        theta_2 = atan2(sqrt(1-OCB(1,1)^2),OCB(1,1)) ;
        theta_3 = atan2(-OCB(1,2),OCB(1,3)) ;
    case 4
        theta_1 = atan2(OCB(3,1),OCB(2,1)) ; 
        theta_2 = atan2(sqrt(1-OCB(1,1)^2),OCB(1,1)) ;
        theta_3 = atan2(OCB(1,3),-OCB(1,2)) ;
    case 5
        theta_1 = atan2(OCB(3,1),OCB(3,3)) ; 
        theta_2 = atan2(-OCB(2,3),sqrt(1-OCB(2,3)^2)) ;
        theta_3 = atan2(OCB(2,1),OCB(2,2)) ;
    case 6
        theta_1 = atan2(-OCB(3,1),OCB(1,1)) ; 
        theta_2 = atan2(OCB(2,1),sqrt(1-OCB(2,1)^2)) ;
        theta_3 = atan2(-OCB(2,3),OCB(2,2)) ;
    case 7
        theta_1 = atan2(OCB(1,2),OCB(3,2)) ; 
        theta_2 = atan2(sqrt(1-OCB(2,2)^2),OCB(2,2)) ;
        theta_3 = atan2(OCB(2,1),-OCB(2,3)) ;
    case 8
        theta_1 = atan2(OCB(3,2),-OCB(1,2)) ; 
        theta_2 = atan2(sqrt(1-OCB(2,2)^2),OCB(2,2)) ;
        theta_3 = atan2(OCB(2,3),OCB(2,1)) ;
    case 9
        theta_1 = atan2(-OCB(1,2),OCB(2,2)) ; 
        theta_2 = atan2(OCB(3,2),sqrt(1-OCB(3,2)^2)) ;
        theta_3 = atan2(-OCB(3,1),OCB(3,3)) ;
    case 10
        theta_1 = atan2(OCB(2,1),OCB(1,1)) ; 
        theta_2 = atan2(-OCB(3,1),sqrt(1-OCB(3,1)^2)) ;
        theta_3 = atan2(OCB(3,2),OCB(3,3)) ;
    case 11
        theta_1 = atan2(OCB(1,3),-OCB(2,3)) ; 
        theta_2 = atan2(sqrt(1-OCB(3,3)^2),OCB(3,3)) ;
        theta_3 = atan2(OCB(3,1),OCB(3,2)) ;
    case 12 
        theta_1 = atan2(OCB(2,3),OCB(1,3)) ; 
        theta_2 = atan2(sqrt(1-OCB(3,3)^2),OCB(3,3)) ;
        theta_3 = atan2(OCB(3,2),-OCB(3,1)) ;
    otherwise
        %fprintf('\n Holy Shit !!! \n');
        disp('Choice is not a valid selection. Returning null angles.')
        theta_1 = 0;
        theta_2 = 0;
        theta_3 = 0;
end
        
% 
% theta_1 = rad2deg(theta_1) ;
% theta_2 = rad2deg(theta_2) ;
% theta_3 = rad2deg(theta_3) ;

end