%Finite element analysis for a circle
%Perimeter of a circle. 
%Assuming uniform meshing.
clc
n = input('Enter the number of elements (meshes) of the circle: ');
rad = input('Enter the radius of the circle: ');

theta = 2*pi/n;

l = 2*rad*sin(theta/2);

perimeter = l*n;

error = abs(perimeter - 2*pi*rad);

fprintf('The perimeter of the circle was found out by using finite element analysis:\n');
disp(perimeter);
fprintf('The error in measuring the perimeter:\n');
disp(error);
    