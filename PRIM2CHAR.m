function [RMAT,LMAT,LAMS] = PRIM2CHAR(qc0,qp1,g,ii)
% Inputs: states at i, i+1 grid pts
% Outputs: right and left eigenvectors, eigenvalues for char variable
% mapping about pt i

rl=qc0(1);
ul=qc0(2)/rl;
El=qc0(3);

rr=qp1(1);
ur=qp1(2)/rr;
Er=qp1(3);

gm1=g-1;

pr=(Er-0.5*rr*ur^2)*gm1;
Hr = Er + pr;
htr = Hr/rr;
% htr = Hr/rr;
% htr=g/gm1*pr/rr + 0.5*ur^2; % MCA THESIS DEFN

pl=(El-0.5*rl*ul^2)*gm1;
Hl = El + pl;
% htl = Hl/rl;
% htl=g/gm1*pl/rl + 0.5*ul^2; % MCA THESIS DEFN
htl=Hl/rl;

ravg=sqrt(rr*rl);
uavg=(sqrt(rr)*ur + sqrt(rl)*ul)/(sqrt(rr)+sqrt(rl));
havg=(sqrt(rr)*htr + sqrt(rl)*htl)/(sqrt(rr)+sqrt(rl));
cavg=sqrt(gm1*(havg-0.5*uavg^2));


if ii==101
  db=1;
end

% Construct matrix of inverse eigenvectors to A = dF/dU
%                          _                                       _
%                         |                                         |
%                         |  uc/(gamma-1)+u^2/2  -c/(gamma-1)-u   1 |
%                         |                                         |
% R^{-1}=(gamma-1)/(2c^2)*|  2(H-u^2)             2u             -2 |
%                         |                                         |
%                         | -uc/(gamma-1)+u^2/2   c/(gamma-1)-u   1 |
%                         |_                                       _|


ic = gm1/(2*cavg^2);

l1 = [uavg*cavg/gm1 + 0.5*uavg^2 2*(havg-uavg^2) -uavg*cavg/gm1 + 0.5*uavg^2]';

l2 = [-cavg/gm1 - uavg  2*uavg  cavg/gm1-uavg]';

l3 = [1 -2 1]';

LMAT=[l1 l2 l3];
LMAT=ic*LMAT;



% Construct matrix classical eigenvectors for A = dF/dU
%      _                    _
%     |                      |
%     |   1      1       1   |
%     |                      |
% R = |  u-c     u      u+c  |
%     |                      |
%     |  H-uc   u^2/2   H+uc |
%     |_                    _|
rc=ravg/(2*cavg);

r1 = [1  uavg-cavg  havg-uavg*cavg]';
r2 = [1 uavg 0.5*uavg^2]';
r3 = [1 uavg+cavg havg+uavg*cavg]';

RMAT = [r1 r2 r3];

% LMAT=inv(RMAT);

%%% ==== READ THIS === %%%
% fprintf(1,'-------- IMPORTANT --------\n')
% fprintf(1,'NEED TO ENSURE EACH EIGENVECTOR AND EIGENVALUE IS LINED UP PROPERLY IN THIS PRIM2CHAR ROUNTINE\n')
% fprintf(1,'-------- IMPORTANT --------\n')
% fprintf(1,'LOOK AT CALEY BOOK OR SOMETHING')
% 3 1D Euler equation eigenvalues

% === NEED TO ENTROPY FIX THE EIGENVALUES?? === %
lam1=uavg-cavg;
lam2=uavg;
lam3=uavg+cavg;


%% Entropy fix
eps=.01;
lamvec=[lam1 lam2 lam3];
mxwv=abs(uavg)+cavg;
mxwv=max(abs(lamvec));
h=eps*mxwv;
del=h;

% if abs(lam1)<del
%   l1fix=(lam1^2+ del^2)/(2*del);
%   fprintf(1,'Point ii = %i  Fixing left acoustic wave\n',ii)
%   fprintf(1,'Old = %f and New = %f\n',lam1,l1fix)
%   lam1=l1fix;
% end

% if abs(lam2)<del
%   %   l2fix=0.5*((lam2^2)/del + del^2);
%   l2fix=(lam2^2+ del^2)/(2*del);
%   fprintf(1,'Point ii=%i   Fixing entropy wave\n',ii)
%   fprintf(1,'Old = %f and New = %f\n',lam2,l2fix)
%   lam2=l2fix;
% end

% if abs(lam3)<del
%   %     l3fix=sign(lam3)*0.5*((lam3^2)/del + del);
%   l3fix=(lam3^2+ del^2)/(2*del);
%   fprintf(1,'Point ii = %i  Fixing right acoustic wave\n',ii)
%   fprintf(1,'Old = %f and New = %f\n',lam3,l3fix)
%   lam3=l3fix;
% end



LAMS=[lam1 lam2 lam3];


% minlam=min(abs(LAMS));
% maxlam=max(abs(LAMS));
% if minlam<eps
% imin=find(abs(LAMS)==minlam);
% imax=find(abs(LAMS)==maxlam);
% LAMS(imin) = 0.05*LAMS(imax(1));


end
