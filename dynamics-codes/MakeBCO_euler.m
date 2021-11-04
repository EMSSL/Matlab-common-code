function [BCO] = MakeBCO_euler(ang1, ang2, ang3, seq)
% Code to make a BCO matrix using numerical values. 
%   INPUTS:
%       ang1, ang2, ang3 = numeric values (rad) of angles
%       seq              = rotation scheme
%               'XYZ' yields vec_B = [RZ*RY*RX] * vec_O
% Christopher D. Yoder
% 
% 
% 2020-10-18, initial revision
% 2021-06-14, rev1, fixed issues with symbolics v floats 

% XYZ -> BCO = MakeBCO_euler('psi', 'phi', 'theta', 'XYZ')

% make matrices
R1 = populate(seq(1), ang1);
if isa(ang1, 'float') == 1
    if mod(ang1, pi/2) == 0
        for i1 = 1:9
            R1(i1) = int8(R1(i1));
        end
    end
end
R2 = populate(seq(2), ang2);
if isa(ang2, 'float') == 1
    if mod(ang2, pi/2) == 0
        for i1 = 1:9
            R2(i1) = int8(R2(i1));
        end
    end
end
R3 = populate(seq(3), ang3);
if isa(ang3, 'float') == 1
    if mod(ang3, pi/2) == 0
        for i1 = 1:9
            R3(i1) = int8(R3(i1));
        end
    end
end

% multiply - euler sequence notation
% EulerSequence(X, Y, Z) -> BCO = RZ*RY*RX
% keyboard
BCO = R3*R2*R1;


end
