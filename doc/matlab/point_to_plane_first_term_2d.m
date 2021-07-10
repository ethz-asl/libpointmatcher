clc;
clear; 
% declare the symbolic variables
x = sym('x','real');
y = sym('y','real');

a = sym('a','real');


T = [x; y;];
R = [cos(a) -sin(a);
     sin(a) cos(a)];

pix = sym('pix','real');
piy = sym('piy','real');

qix = sym('qix','real');
qiy = sym('qiy','real');
 
Pi = [pix;piy];
Qi = [qix;qiy];

nix = sym('nix','real');
niy = sym('niy','real');

Ni = [nix;niy;];

G = dot((R * Pi + T - Qi), Ni);    %Ni should correspond to Qi
dG_dx = jacobian(G,x);
dJ_dx = 2 * dG_dx' * G;

dG_dy = jacobian(G,y);
dJ_dy = 2 * dG_dy' * G;


%% 1,2
d2J_dx2 = jacobian(dJ_dx,x)
d2J_dy2 = jacobian(dJ_dy,y)
%% 4,5,6 , 34,35,36
d2J_dydx = jacobian(dJ_dx,y)
d2J_dxdy = jacobian(dJ_dy,x)




%% 7
dG_da = jacobian(G,a);
dJ_da = 2 * dG_da' * G;
d2J_da2 = jacobian(dJ_da,a)


%% (10,11) , (12,13) , (14,15)
d2J_dxda = jacobian(dJ_da,x)
d2J_dadx = jacobian(dJ_dx,a)

d2J_dyda = jacobian(dJ_da,y)
d2J_dady = jacobian(dJ_dy,a)





































