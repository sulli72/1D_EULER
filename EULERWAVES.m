function [lam1,lam2,lam3] = EULERWAVES(q,g)


[r,u,p] = Q2PRIM(q,g);

c=sqrt(g*p/r);

lam1=u+c;
lam2=u;
lam3=u-c;

% eps=0.05*lam1;
% lam3=max(abs(lam3),abs(eps));

end