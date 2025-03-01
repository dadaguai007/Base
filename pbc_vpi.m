%PBC
function [Ex_out, Ey_out]= pbc_vpi(Ex,Ey,theta)
%Not tested ， the function is done
if nargin < 2
    theta = 0;  % Default rotation angle
end

% Check and reshape E if necessary
if size(Ex, 2) == 1
    Ex = [Ex, zeros(size(Ex))];  % Assume Y component is 0
end
% Check and reshape E if necessary
if size(Ey, 2) == 1
    Ey = [zeros(size(Ey)),Ey];  % Assume Y component is 0
end
theta1=theat+pi/2;
J1=[cos(theta).^2,cos(theta)*sin(theta);cos(theta)*sin(theta),sin(theta).^2];
J2=[cos(theta1).^2,cos(theta1)*sin(theta1);cos(theta1)*sin(theta1),sin(theta1).^2];
Eout1=Ex*J1;
Eout2=Ey*J2;
E=Eout1+Eout2;
Ex_out = E(:, 1);
Ey_out = E(:, 2);
end