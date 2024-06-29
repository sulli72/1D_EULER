function FLX = Q2FLUX(Q,g)
% This function takes the solution Q at each grid point and computes the 
% flux function F(Q) at each grid point.

gm1=g-1;

[r,u,p] = Q2PRIM(Q,g);

E= p/(gm1) + 0.5*r.*u.^2; % Compute total energy

f1 = r.*u; % <--- mass flux
f2 = r.*u.^2+p; % <--- momentum flux
f3 = u.*(E+p); % <--- energy flux 
% (E+p is equivalent to denisty * total enthalpy/ unit mass 
% i.e E+p = H = r*h = r*(e+p/r) 

FLX = [f1 f2 f3];

end