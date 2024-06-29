function [VPLS,VMNS] = CHARFLUX(fpls,fmns,Q,g)

nx=size(Q,1);
neqs=size(Q,2);
VPLS = zeros(nx,neqs);
VMNS = zeros(nx,neqs);


for ii=1:nx

  % Get left and right states for Roe avg
  % --- same eigenvectors for +/- flux (bc frozen at x+1/2 pts)
  qc0=Q(ii,:);
  if ii==nx
    qp1=qc0;
  else
    qp1=Q(ii+1,:);
  end

  myfp = fpls(ii,:);
  myfm = fmns(ii,:);


   rr=qp1(1);   rl=qc0(1);
   ur=qp1(2);   ul=qc0(2);
   Er=qp1(3);   El=qc0(3);

   gm1=g-1;

   pr = gm1*(Er-0.5*rr*ur*ur);
   Hr = Er + pr;
   htr = Hr/rr;

   pl = gm1*(El-0.5*rl*ul*ul);
   Hl = El + pl;
   htl = Hl/rl;

   ravg=sqrt(rr*rl);

   uavg=(sqrt(rr)*ur + sqrt(rl)*ul)/(sqrt(rr)+sqrt(rl));

   havg=(sqrt(rr)*htr + sqrt(rl)*htl)/(sqrt(rr)+sqrt(rl));

   cavg=sqrt((g-1)*(havg-0.5*uavg^2));

  % Right characteristic (eigen) vectors
  cr2=ravg/(2*cavg); %coeffs for r2,r3
  cr3=-cr2;
  
  r1=[1 uavg 0.5*uavg^2]';

%   r2=[1 uavg+cavg havg+cavg*uavg]';
  r2 = [1 uavg+cavg uavg^2/2 + cavg^2*gm1 + uavg*cavg]';
  r2=cr2*r2;

%   r3=[1 uavg-cavg havg-cavg*uavg]';
  r3 = [1 uavg+cavg uavg^2/2 + cavg^2*gm1 - uavg*cavg]';
  r3=cr3*r3;
  
  RMAT=[r1 r2 r3];

  ic = (gm1)/(ravg*cavg);
  rbc=ravg/cavg;
  cuugm1=cavg*uavg^2/gm1;
  uuh=uavg*uavg*0.5;
  cgm1=cavg/gm1;
  ccgm1=cavg*cavg/gm1;

  l1 = [rbc*(ccgm1-uuh) uuh-cuugm1 -uuh-cuugm1]';
  l1 = ic*l1;

  l2 = [rbc*uavg cgm1-uavg cgm1+uavg]';
  l2 = ic*l2;

  l3 = [-rbc 1 -1]';
  l3 = ic*l3;

  LMAT = [l1 l2 l3];
%   LMAT=LMAT';

%   LMAT = inv(RMAT);

%   RINV=inv(RMAT);
  VPLS(ii,:) = LMAT*myfp';
  VMNS(ii,:) = LMAT*myfm';

  db=1;




end

end