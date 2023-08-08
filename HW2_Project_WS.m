%%
%WorkSpace
% Inserting D-H convention parameters
a1 = 0; alpha1 = -pi/2; d1 = 30;
a2 = 0; alpha2 = -pi/2; d2 = 0;
a3 =-0.045; alpha3 = pi/2; t3 = 0;

% Inserting joint limits for Arms
t1_min = -pi/4; t1_max = pi/4;
t2_min = -pi/4; t2_max = pi/4;
d3_min = 0; d3_max = 3;

% Monte Carlo method
% sampling size
N = 30000;
t1 = t1_min + (t1_max-t1_min)*rand(N,1);
t2 = t2_min + (t2_max-t2_min)*rand(N,1);
d3 = d3_min + (d3_max-d3_min)*rand(N,1);




for i = 1:N
A1 = TransMat(a1,alpha1,d1,t1(i));
A2 = TransMat(a2,alpha2,d2,t2(i));
A3 = TransMat(a3,alpha3,d3(i),t3);

T = A1*A2*A3;
X=T(1,4);
Y=T(2,4);
Z=T(3,4);
plot3(X,Y,Z,'.')
hold on;
end

function [ T ] = TransMat( a,b,c,d )

T = [ cos(d) -sin(d)*cos(b) sin(d)*sin(b) a*cos(d);
    sin(d) cos(d)*cos(b) -cos(d)*sin(b) a*sin(d);
0 sin(b) cos(b) c;
0 0 0 1];
end

