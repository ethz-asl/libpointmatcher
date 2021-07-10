clc;
clear; 
% declare the symbolic variables
x = sym('x','real');
y = sym('y','real');
z = sym('z','real');

a = sym('a','real');
b = sym('b','real');
c = sym('c','real');


T = [x; y; z;];
R = [cos(a)*cos(b) cos(a)*sin(b)*sin(c)-sin(a)*cos(c) cos(a)*sin(b)*cos(c)+sin(a)*sin(c);
     sin(a)*cos(b) sin(a)*sin(b)*sin(c)+cos(a)*cos(c) sin(a)*sin(b)*cos(c)-cos(a)*sin(c);
     -sin(b)       cos(b)*sin(c)                      cos(b)*cos(c)                      ];

pix = sym('pix','real');
piy = sym('piy','real');
piz = sym('piz','real'); 

qix = sym('qix','real');
qiy = sym('qiy','real');
qiz = sym('qiz','real');
 
Pi = [pix;piy;piz];
Qi = [qix;qiy;qiz];

nix = sym('nix','real');
niy = sym('niy','real');
niz = sym('niz','real');

Ni = [nix;niy;niz];

G = dot((R * Pi + T - Qi), Ni);    %Ni should correspond to Qi
dG_dx = jacobian(G,x);
dJ_dx = 2 * dG_dx' * G;

dG_dy = jacobian(G,y);
dJ_dy = 2 * dG_dy' * G;

dG_dz = jacobian(G,z);
dJ_dz = 2 * dG_dz' * G;

%% 1,2,3
d2J_dx2 = jacobian(dJ_dx,x)
d2J_dy2 = jacobian(dJ_dy,y)
d2J_dz2 = jacobian(dJ_dz,z)
%% 4,5,6 , 34,35,36
d2J_dydx = jacobian(dJ_dx,y)
d2J_dxdy = jacobian(dJ_dy,x)

d2J_dzdx = jacobian(dJ_dx,z)
d2J_dxdz = jacobian(dJ_dz,x)

d2J_dydz = jacobian(dJ_dz,y)
d2J_dzdy = jacobian(dJ_dy,z)

%% 7
dG_da = jacobian(G,a);
dJ_da = 2 * dG_da' * G;
d2J_da2 = jacobian(dJ_da,a)
%% 8
dG_db = jacobian(G,b);
dJ_db = 2 * dG_db' * G;
d2J_db2 = jacobian(dJ_db,b)
%% 9
dG_dc = jacobian(G,c);
dJ_dc = 2 * dG_dc' * G;
d2J_dc2 = jacobian(dJ_dc,c)
%% (10,11) , (12,13) , (14,15)
d2J_dxda = jacobian(dJ_da,x)
d2J_dadx = jacobian(dJ_dx,a)

d2J_dyda = jacobian(dJ_da,y)
d2J_dady = jacobian(dJ_dy,a)

d2J_dzda = jacobian(dJ_da,z)
d2J_dadz = jacobian(dJ_dz,a)


%% (16,17), (18,19), (20,21)
d2J_dxdb = jacobian(dJ_db,x)
d2J_dbdx = jacobian(dJ_dx,b)

d2J_dydb = jacobian(dJ_db,y)
d2J_dbdy = jacobian(dJ_dy,b)

d2J_dzdb = jacobian(dJ_db,z)
d2J_dbdz = jacobian(dJ_dz,b)

%% (22,23) ,(24,25) , (26,27)
d2J_dxdc = jacobian(dJ_dc,x)
d2J_dcdx = jacobian(dJ_dx,c)

d2J_dydc = jacobian(dJ_dc,y)
d2J_dcdy = jacobian(dJ_dy,c)

d2J_dzdc = jacobian(dJ_dc,z)
d2J_dcdz = jacobian(dJ_dz,c)

%% 28,29 , 30,31,  32,33

d2J_dadb = jacobian(dJ_db,a)
d2J_dbda = jacobian(dJ_da,b)

d2J_dbdc = jacobian(dJ_dc,b)
d2J_dcdb = jacobian(dJ_db,c)

d2J_dcda = jacobian(dJ_da,c)
d2J_dadc = jacobian(dJ_dc,a)

































