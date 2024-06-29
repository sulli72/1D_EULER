function [r,u,p] = CONS2PRIM(R,RU,E,g)
% Subroutine to convert conservative variables into primitive form

gm1=g-1;
nx=length(R);
r=zeros(nx,1);
u=zeros(nx,1);
p=zeros(nx,1);

r=R;
u=RU./R;
p=(E-0.5*r.*u.^2)*gm1;

end