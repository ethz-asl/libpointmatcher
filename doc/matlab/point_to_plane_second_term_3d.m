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
%% small trick
dG_dT = jacobian(G,T);
dJ_dT = 2 * dG_dT' * G;
d2J_dT2 = jacobian(dJ_dT,T);

%%
dG_dx = jacobian(G,x);
dG_dy = jacobian(G,y);
dG_dz = jacobian(G,z);
dG_da = jacobian(G,a);
dG_db = jacobian(G,b);
dG_dc = jacobian(G,c);

dJ_dx = 2 * dG_dx' * G;
dJ_dy = 2 * dG_dy' * G;
dJ_dz = 2 * dG_dz' * G;
dJ_da = 2 * dG_da' * G;
dJ_db = 2 * dG_db' * G;
dJ_dc = 2 * dG_dc' * G;

%% Pi
d2J_dpix_dx = jacobian(dJ_dx, pix)
d2J_dpix_dy = jacobian(dJ_dy, pix)
d2J_dpix_dz = jacobian(dJ_dz, pix)
d2J_dpix_da = jacobian(dJ_da, pix)
d2J_dpix_db = jacobian(dJ_db, pix)
d2J_dpix_dc = jacobian(dJ_dc, pix)

d2J_dpiy_dx = jacobian(dJ_dx, piy)
d2J_dpiy_dy = jacobian(dJ_dy, piy)
d2J_dpiy_dz = jacobian(dJ_dz, piy)
d2J_dpiy_da = jacobian(dJ_da, piy)
d2J_dpiy_db = jacobian(dJ_db, piy)
d2J_dpiy_dc = jacobian(dJ_dc, piy)

d2J_dpiz_dx = jacobian(dJ_dx, piz)
d2J_dpiz_dy = jacobian(dJ_dy, piz)
d2J_dpiz_dz = jacobian(dJ_dz, piz)
d2J_dpiz_da = jacobian(dJ_da, piz)
d2J_dpiz_db = jacobian(dJ_db, piz)
d2J_dpiz_dc = jacobian(dJ_dc, piz)

%% Qi
d2J_dqix_dx = jacobian(dJ_dx, qix)
d2J_dqix_dy = jacobian(dJ_dy, qix)
d2J_dqix_dz = jacobian(dJ_dz, qix)
d2J_dqix_da = jacobian(dJ_da, qix)
d2J_dqix_db = jacobian(dJ_db, qix)
d2J_dqix_dc = jacobian(dJ_dc, qix)

d2J_dqiy_dx = jacobian(dJ_dx, qiy)
d2J_dqiy_dy = jacobian(dJ_dy, qiy)
d2J_dqiy_dz = jacobian(dJ_dz, qiy)
d2J_dqiy_da = jacobian(dJ_da, qiy)
d2J_dqiy_db = jacobian(dJ_db, qiy)
d2J_dqiy_dc = jacobian(dJ_dc, qiy)

d2J_dqiz_dx = jacobian(dJ_dx, qiz)
d2J_dqiz_dy = jacobian(dJ_dy, qiz)
d2J_dqiz_dz = jacobian(dJ_dz, qiz)
d2J_dqiz_da = jacobian(dJ_da, qiz)
d2J_dqiz_db = jacobian(dJ_db, qiz)
d2J_dqiz_dc = jacobian(dJ_dc, qiz)













