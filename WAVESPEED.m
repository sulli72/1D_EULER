function speed = WAVESPEED(Q,g)
% This function determins the maximum wavespeed in the 1D Euler equations
% which is given by taking max(abs(u)+c) over all grid points

% State variables
[r,u,p]=Q2PRIM(Q,g);
nx=length(r);

% Speed of sound vector
csq=zeros(nx,1);
csq = g*p./r;
c = csq.^0.5;

wave = abs(u) + c;
speed=max(abs(wave));

end