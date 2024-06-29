function [FPLS,FMNS] = CHAR2PRIMFLUX(vpls,vmns,Q,g)

nx=size(Q,1);
neqs=size(Q,2);
FPLS = zeros(nx,neqs);
FMNS = zeros(nx,neqs);


for ii=1:nx

  % Get left and right states for Roe avg
  % --- same eigenvectors for +/- flux (bc frozen at x+1/2 pts)
  qc0=Q(ii,:);
  if ii==nx
    qp1=qc0;
  else
    qp1=Q(ii+1,:);
  end

  myfp = vpls(ii,:);
  myfm = vmns(ii,:);


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
  r3 = [1 uavg+cavg uavg^2/2 + cavg^2*gm1-uavg*cavg]';
  r3=cr3*r3;
  
  RMAT=[r1 r2 r3];

%   RINV=inv(RMAT);
  FPLS(ii,:) = RMAT*myfp';
  FMNS(ii,:) = RMAT*myfm';

  db=1;








end

end