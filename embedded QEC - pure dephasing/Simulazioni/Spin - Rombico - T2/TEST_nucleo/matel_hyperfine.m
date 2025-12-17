clear all
close all
I = 3/2;
S = 1/2;
dim = (2*S+1)*(2*I+1);

mu_b  = 0.057883826350000e+00;
g = 2;

A = [400 400 900]; % MHz, VO
p = 66; % MHz, 
g = 2;
muB = 14; % MHz/mT
B = 200;  %mT
theta = pi/2;

Ix = stev(I,1, 1);
Iy = stev(I,1, -1);
Iz = stev(I,1, 0);

Sx = stev(S,1, 1);
Sy = stev(S,1, -1);
Sz = stev(S,1, 0);

k = 1;

Ham = A(1)*kron(Sx,Ix) + A(2)*kron(Sy,Iy) + A(3)*kron(Sz,Iz) + ...
      p*kron(eye(2*S+1),Iz^2) + ...
      g*muB*B*( kron(Sx,eye(2*I+1)) * sin(theta) + kron(Sz,eye(2*I+1)) * cos(theta) );

[V,ev] = eig(Ham, 'vector');

Ix_r = V'* kron(eye(2*S+1),Ix) * V;
Iy_r = V'* kron(eye(2*S+1),Iy) * V;
Iz_r = V'* kron(eye(2*S+1),Iz) * V;
Sx_r = V'* kron(Sx,eye(2*I+1)) * V;
Sy_r = V'* kron(Sy,eye(2*I+1)) * V;
Sz_r = V'* kron(Sz,eye(2*I+1)) * V;

I_par = Ix_r * sin(theta) + Iz_r * cos(theta) ;
S_par = Sx_r * sin(theta) + Sz_r * cos(theta) ;
[[1:dim].' diag(I_par) diag(S_par)]

abs(g*mu_b*Sx_r + g*mu_b/1880 * Ix_r) + abs(g*mu_b*Sz_r + g*mu_b/1880 * Iz_r)
