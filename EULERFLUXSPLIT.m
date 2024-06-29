function [FPLS,FMNS] = EULERFLUXSPLIT(FLXVEC,Q,g)

% State variables
gm1=g-1;
[r,u,p]=Q2PRIM(Q,g);
E = p/gm1 + 0.5*r.*u.*u;
nx=length(r);
neqs=size(Q,2);

% Speed of sound vector
c=zeros(nx,1);
c = g*p./r;


lam1 = u + c;
lam2 = u;
lam3 = u - c;

aupc=abs(u)+c;
alpha=max(aupc);
ialph=1/alpha;


%%% GLOBAL LAX FRIEDRICHS SPLITTING
FPLS = 1/2*(FLXVEC + alpha*Q);
FMNS = 1/2*(FLXVEC - alpha*Q);


db=1;
end