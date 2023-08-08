clc
close all
clear
syms theta theta1 alpha d a theta2 theta3 L1 L2 L3;



%% 
alpha = sym(-pi/2);
theta = theta1; 
d=L1+L3;
a=0;
T1=[cos(theta), -sin(theta)*cos(alpha), sin(theta)*sin(alpha), a*cos(theta);
    sin(theta), cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta);
    0         , sin(alpha)           , cos(alpha)            ,            d;
    0         , 0                    , 0                     ,            1];

alpha = sym(-pi/2);
theta = theta2; 
d=0;
a=0;
T2=[cos(theta), -sin(theta)*cos(alpha), sin(theta)*sin(alpha), a*cos(theta);
    sin(theta), cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta);
    0         , sin(alpha)           , cos(alpha)            ,            d;
    0         , 0                    , 0                     ,            1];
                
alpha = 0;
theta = 0; 
d=d;
a=0;
T3=[cos(theta), -sin(theta)*cos(alpha), sin(theta)*sin(alpha), a*cos(theta);
    sin(theta), cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta);
    0         , sin(alpha)           , cos(alpha)            ,            d;
    0         , 0                    , 0                     ,            1];   
                
T=simplify(combine(T1*T2*T3,'sincos'))


%%
syms px py pz sx sy sz sox soy soz ;
st=sin(theta);
ct=cos(theta);
vt=1-cos(theta);
d=sox*sx+soy*sy+soz*sz;

px=d*sx-sox*(sx^2-1)*vt-soy*(sx*sy*vt-sz*st)-soz*(sx*sz*vt+sy*st);
py=d*sy-sox*(sy*sx*vt+sz*st)-soy*(sy^2-1)*vt-soz*(sy*sz*vt-sx*st);
pz=d*sz-sox*(sz*sx*vt-sy*st)-soy*(sz*sy*vt+sx*st)-soz*(sz^2-1)*vt;

A=[sx^2*vt+ct, sx*sy*vt-sz*st, sx*sz*vt+sy*st, px;
   sy*sx*vt+sz*st, sy^2*vt+ct, sz*sy*vt-sx*st, py;
   sx*sz*vt-sy*st, sy*sz*vt+sx*st, sz^2*vt+ct,pz;
   0             , 0             , 0         , 1];


%%
sx=0;
sy=0;
sz=1;
sox=0;
soy=0;
soz=0;
st=sin(theta1);
ct=cos(theta1);
vt=1-cos(theta1);
d=sox*sx+soy*sy+soz*sz;

px=d*sx-sox*(sx^2-1)*vt-soy*(sx*sy*vt-sz*st)-soz*(sx*sz*vt+sy*st);
py=d*sy-sox*(sy*sx*vt+sz*st)-soy*(sy^2-1)*vt-soz*(sy*sz*vt-sx*st);
pz=d*sz-sox*(sz*sx*vt-sy*st)-soy*(sz*sy*vt+sx*st)-soz*(sz^2-1)*vt;
A1=[sx^2*vt+ct, sx*sy*vt-sz*st, sx*sz*vt+sy*st, px;
   sy*sx*vt+sz*st, sy^2*vt+ct, sz*sy*vt-sx*st, py;
   sx*sz*vt-sy*st, sy*sz*vt+sx*st, sz^2*vt+ct,pz;
   0             , 0             , 0         , 1];

%%
sx=0;
sy=1;
sz=0;
sox=0;
soy=0;
soz=L1+L3;
st=sin(theta2);
ct=cos(theta2);
vt=1-cos(theta2);
d=sox*sx+soy*sy+soz*sz;

px=d*sx-sox*(sx^2-1)*vt-soy*(sx*sy*vt-sz*st)-soz*(sx*sz*vt+sy*st);
py=d*sy-sox*(sy*sx*vt+sz*st)-soy*(sy^2-1)*vt-soz*(sy*sz*vt-sx*st);
pz=d*sz-sox*(sz*sx*vt-sy*st)-soy*(sz*sy*vt+sx*st)-soz*(sz^2-1)*vt;
A2=[sx^2*vt+ct, sx*sy*vt-sz*st, sx*sz*vt+sy*st, px;
   sy*sx*vt+sz*st, sy^2*vt+ct, sz*sy*vt-sx*st, py;
   sx*sz*vt-sy*st, sy*sz*vt+sx*st, sz^2*vt+ct,pz;
   0             , 0             , 0         , 1];

%%
sx=0;
sy=1;
sz=0;
sox=0;
soy=0;
soz=L1+L3;
theta3=0;
st=sin(theta3);
ct=cos(theta3);
vt=1-cos(theta3);
d=sox*sx+soy*sy+soz*sz;

px=d*sx-sox*(sx^2-1)*vt-soy*(sx*sy*vt-sz*st)-soz*(sx*sz*vt+sy*st);
py=d*sy-sox*(sy*sx*vt+sz*st)-soy*(sy^2-1)*vt-soz*(sy*sz*vt-sx*st);
pz=d*sz-sox*(sz*sx*vt-sy*st)-soy*(sz*sy*vt+sx*st)-soz*(sz^2-1)*vt;
A3=[sx^2*vt+ct, sx*sy*vt-sz*st, sx*sz*vt+sy*st, px;
   sy*sx*vt+sz*st, sy^2*vt+ct, sz*sy*vt-sx*st, py;
   sx*sz*vt-sy*st, sy*sz*vt+sx*st, sz^2*vt+ct,pz;
   0             , 0             , 0         , 1];


AM=A1*A2*A3*[1, 0 , 0, 0; 0, -1, 0, 0; 0, 0, -1, L1+L3-d; 0, 0, 0, 1];

A=simplify(combine(AM,'sincos'))




%%
J=inv(T2)*inv(T1)
L=simplify(combine(J,'sincos'))
T3

W=simplify(combine(inv(T2),'sincos'))




