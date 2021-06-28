function [BCO] = MakeBCO_quat(q0, q1, q2, q3)

% Written by Sathvick Divi for EMSSL on 05/19/2016. Modified by Christopher D. Yoder for EMSSL
% on 05/21/2016. 
%
% Modified on 11/11/2016 by Christopher D. Yoder and Jacob D. Reedy to
% correct errors in the OCB computation.

% also a built-in matlab to normalize quaterions called 'quatnomrlaize'
% but requires the aerospace toolbox from matlab

% check symbolic or numeric
% keyboard
if isa(q0, 'symfun') || isa(q0, 'sym')
% if strcmp(class(q0),'sym') == 1         % check for symbolic variables
    BCO = sym('BCO',[3,3]);     
    q0t = q0;
    q1t = q1;
    q2t = q2;
    q3t = q3;
elseif isnumeric(q0) == 1               % check for numeric variables
    BCO = zeros(3, 3);
    [q0t, q1t, q2t, q3t] = quat_normalize(q0, q1, q2, q3);
else
    warning('q0 is neither symbolic or numeric. Attempting to execute anyway.');
end

% assemble the BCO matrix
% BCO(1,1) = q0t^2 + q1t^2 + q2t^2 + q3t^2 ;
BCO(1,1) = q0t^2 + q1t^2 - q2t^2 - q3t^2 ;
BCO(1,2) = 2*(q1t*q2t + q0t*q3t) ;
BCO(1,3) = 2*(q1t*q3t - q0t*q2t) ;
BCO(2,1) = 2*(q1t*q2t - q0t*q3t) ;
% BCO(2,2) = (q0t^2 - q1t^2) + (q2t^2 + q3t^2) ;
BCO(2, 2) = q0t^2 - q1t^2 + q2t^2 - q3t^2 ;
BCO(2,3) = 2*(q2t*q3t + q0t*q1t);
BCO(3,1) = 2*(q1t*q3t + q0t*q2t);
BCO(3,2) = 2*(q2t*q3t - q0t*q1t);
BCO(3,3) = q0t^2 - q1t^2 - q2t^2 + q3t^2 ;

end