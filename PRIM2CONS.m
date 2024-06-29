function [R,RU,E] = PRIM2CONS(r,u,p,g)
%%% This function turns the primitive variables [r,u,p] into their
%%% conservative form [r,ru,E] 
gm1=g-1;
R=r;
RU=r.*u;
E= p/(gm1) + 0.5*r.*u.^2;
end