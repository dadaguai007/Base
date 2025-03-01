A=[1,0;0,0];
B=[0,0;0,1];
E=ones(10,1);
E2=zeros(10,1);
E=[E2,E];
E1=[0,1];
theta=-pi/4;
J=[cos(theta).^2,cos(theta)*sin(theta);cos(theta)*sin(theta),sin(theta).^2];

% PBC

E1=[Ex1,Ey1];
E2=[Ex2,Ey2];
theta=0;
theta1=theat+pi/2;
J1=[cos(theta).^2,cos(theta)*sin(theta);cos(theta)*sin(theta),sin(theta).^2];
J2=[cos(theta1).^2,cos(theta1)*sin(theta1);cos(theta1)*sin(theta1),sin(theta1).^2];
Eout1=E1*J1;
Eout2=E2*J2;
E=Eout1+Eout2;
Ex = E(:, 1);
Ey = E(:, 2);