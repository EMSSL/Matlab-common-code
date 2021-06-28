function [q0n, q1n, q2n, q3n] = quat_normalize(q0, q1, q2, q3)

% error handling
if isnumeric(q0) ~= 1
    warning('q0 value is non-numeric. Attempting to normalize...');
elseif isnumeric(q1) ~= 1
    warning('q1 value is non-numeric. Attempting to normalize...');
elseif isnumeric(q2) ~= 1
    warning('q2 value is non-numeric. Attempting to normalize...');
elseif isnumeric(q3) ~= 1
    warning('q3 value is non-numeric. Attempting to normalize...');
end

% compute the magnitude of the non-normalized quaternion vector
mag = sqrt(q0^2 + q1^2 + q2^2 + q3^2);

% compute the new values of the normalized quaternion vector
q0n = q0/mag;
q1n = q1/mag;
q2n = q2/mag;
q3n = q3/mag;

end