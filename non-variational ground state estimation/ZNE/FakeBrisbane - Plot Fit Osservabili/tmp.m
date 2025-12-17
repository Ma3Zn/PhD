clear
clc

Z = [1 0; 0 -1];
X = [0 1; 1  0];

Rzz = expm(-1i*pi/2*kron(Z,Z));
Rxx = expm(-1i*pi/2*kron(X,X));

I = eye(4);

A = kron(Rzz, eye(4));
B = kron(eye(2), kron(Rxx, eye(2)));
C = kron(eye(4), Rzz);

norm((A+C)*B - B*(A+C))