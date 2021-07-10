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


%% small trick
dG_dT = jacobian(G,T);
dJ_dT = 2 * dG_dT' * G;
d2J_dT2 = jacobian(dJ_dT,T);

%%
dG_dx = jacobian(G,x);
dG_dy = jacobian(G,y);
dG_da = jacobian(G,a);

dJ_dx = 2 * dG_dx' * G;
dJ_dy = 2 * dG_dy' * G;
dJ_da = 2 * dG_da' * G;

%% Pi
d2J_dpix_dx = jacobian(dJ_dx, pix)
d2J_dpix_dy = jacobian(dJ_dy, pix)
d2J_dpix_da = jacobian(dJ_da, pix)

d2J_dpiy_dx = jacobian(dJ_dx, piy)
d2J_dpiy_dy = jacobian(dJ_dy, piy)
d2J_dpiy_da = jacobian(dJ_da, piy)


%% Qi
d2J_dqix_dx = jacobian(dJ_dx, qix)
d2J_dqix_dy = jacobian(dJ_dy, qix)
d2J_dqix_da = jacobian(dJ_da, qix)

d2J_dqiy_dx = jacobian(dJ_dx, qiy)
d2J_dqiy_dy = jacobian(dJ_dy, qiy)
d2J_dqiy_da = jacobian(dJ_da, qiy)















