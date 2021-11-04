% testing rotations, etc. 
% christopher d. yoder
% 06/27/2021

clc
close all
clear 
fclose('all');
% UsePackage('remove', 'DynPackage');
r2d = 180/pi;
d2r = 1/r2d;

% possible rotations
% XYZ, YXZ, ZXY, XZY, YZX, ZYX
% XYX, YXY, ZXZ, XZX, YZY, ZYZ
syms phi psi theta
% seq = 'XYZ'
% seq = 'YXZ'
% seq = 'ZXY'
% seq = 'XZY'
% seq = 'YZX'
seq = 'ZYX'



% math things 
for i1 = 1:3
    switch seq(i1)
        case 'X'
            angi = theta;
        case 'Y'
            angi = phi;
        case 'Z'
            angi = psi;
    end
    eval(sprintf('%s = %s;', sprintf('ang%d', i1), angi));
end
BCO = MakeBCO_euler(ang1, ang2, ang3, seq) % make the euler matrix 



% convert to latex notation
row1 = sprintf('\t%s & %s & %s', BCO(1, 1), BCO(1, 2), BCO(1, 3));
row2 = sprintf('\t%s & %s & %s', BCO(2, 1), BCO(2, 2), BCO(2, 3));
row3 = sprintf('\t%s & %s & %s', BCO(3, 1), BCO(3, 2), BCO(3, 3));

% scan, find
row1ss = strrep(strrep(row1, '*', ''), '(', '(\');
row2ss = strrep(strrep(row2, '*', ''), '(', '(\');
row3ss = strrep(strrep(row3, '*', ''), '(', '(\');

% print to cmd line
fprintf('\\begin{equation*}\n');
fprintf('\t\\RotateMat{B}{O} = \\Rotate{\\%s}{%s}\\Rotate{\\%s}{%s}\\Rotate{\\%s}{%s} \\\\ \n', ang3, seq(3), ang2, seq(2), ang1, seq(1));
fprintf('\\end{equation*}\n');
fprintf('\\begin{equation*}\n');
fprintf('\t\\begin{bmatrix}\n');
fprintf('\t%s \\\\ \n', row1ss);
fprintf('\t%s \\\\ \n', row2ss);
fprintf('\t%s \n', row3ss);
fprintf('\t\\end{bmatrix}\n');
fprintf('\\end{equation*}\n');