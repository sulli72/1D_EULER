function [FIP,FIM] = ROEFLUX(Q,g)
% Get Roe Flux at the i+1/2, i-1/2 points
% Start with first order Roe
% gives at point i+1/2, Ul = U_i , Ur = U_i+1
%
% Goal is to get conservative Finite Difference approach:
% u_i^n+1 = u_i^n - k/h(F_i+1/2 - F_i-1/2)
%
% where F_i+1/2 = F(Ul_i+1/2, Ur_i+1/2)

%                              Cell                                 
%                 ------------------------------
%                 |                            |
%                 |                            |
% F(u)_i-1/2 ---> | i-1/2       i        i+1/2 | F(u)_i+1/2 --->   
%                 |                            |
%                 |                            |
%                 ------------------------------

% OK to use 1st order Roe scheme, even tho this is a finite difference code
% this is b/c node and cell average values are O(dx^2) approximations of
% each other, therefor the error of the Roe scheme itself , O(dx^1),
% outweighs any error incurred using a FV approach in a FD setting


nvar=3;
nx=size(Q,1);
FIP=zeros(nvar,nx); % Flux at i+1/2 points
FIM=zeros(nvar,nx); % Flux at i-1/2 points

%Get prim variables for ease
[r,u,p]=Q2PRIM(Q,g);
gm1=g-1;
E = Q(:,3);
% Compute total enthalpy - used heavily in Roe scheme
H=E+p;
ht=H./r;


%%% Compute Roe Avg States from i=1/2 to i=xf-1/2 points (interfaces)

for ii=1:nx-1
  rr=r(ii+1);   rl=r(ii);
  ur=u(ii+1);   ul=u(ii);
  pr=p(ii+1);   pl=p(ii);
  htr=ht(ii+1); htl=ht(ii);

  ravg(ii)=sqrt(rr*rl);

  uavg(ii)=(sqrt(rr)*ur + sqrt(rl)*ul)/(sqrt(rr)+sqrt(rl));

  havg(ii)=(sqrt(rr)*htr + sqrt(rl)*htl)/(sqrt(rr)+sqrt(rl));

  cavg(ii)=sqrt( gm1*( havg(ii)-0.5*uavg(ii)^2 ) );

  % Differences in primitive variables across cell interface
  drho(ii)=rr-rl;
  dp(ii)=pr-pl;
  du(ii)=ur-ul;


end
db=1;

% Begin Approximate Riemann Solver

for ii=1:nx-1 % index matches length of Roe avg state vectors
  % Get local average values and state variable jumps
  myu=uavg(ii);
  myc=cavg(ii);
  myr=ravg(ii);
  myh=havg(ii);

  mydu=du(ii); mydp=dp(ii); mydrho=drho(ii);

  % Characterstic wave speeds (eigenvalues of Roe averaged Jacobian)
  lam1=myu;
  lam2=myu+myc;
  lam3=myu-myc;
 
  %%% Entropy fix for transonic case -- seems buggy
  eps=.01;
  lams=[lam1 lam2 lam3];
  mxwv=max(abs(lams));  
  del=eps*mxwv;

  if abs(lam2)<del
%   l2fix=0.5*((lam2^2)/del + del^2);
    l2fix=(lam2^2+ del^2)/(2*del);
  lam2=l2fix;
  fprintf(1,'Point ii=%i   Fixing right acoustic wave\n',ii)
  end

  if abs(lam3)<del
%     l3fix=sign(lam3)*0.5*((lam3^2)/del + del);
    l3fix=(lam3^2+ del^2)/(2*del);
    fprintf(1,'Point ii = %i  Fixing left acoustic wave\n',ii)
    fprintf(1,'Old = %f and New = %f\n',lam3,l3fix)
    lam3=l3fix;
  end


  lams=[lam1 lam2 lam3];

  % Characteristic wave strengths
  dv1=mydrho - mydp/myc^2;   % Entropy wave
  dv2=mydu + mydp/(myr*myc); % Right acoustic wave
  dv3=mydu - mydp/(myr*myc); % Left acoustic wave
  WVMAG = [dv1 dv2 dv3]';

  % Compute right characteristic (eigen) vectors of Roe flux Jacobian 
  % R = span([r1 r2 r3])

  cr2=myr/(2*myc); %coeffs for r2,r3
  cr3=-cr2;

  % Construct matrix of eigenvectors
  r1=[1 myu 0.5*myu^2]';
  r2=[1 myu+myc myh+myc*myu]';
  r2=cr2*r2;
  r3=[1 myu-myc myh-myc*myu]';
  r3=cr3*r3;
  RMAT=[r1 r2 r3];


  % Compute fluxes via combination of left state (ul,rl,pl) and appropriate
  % linear combination of eigenvectors for flux jacobian
  % THIS IS ALLOWED BC ROE SCHEME LINEARIZES THE HYPERBOLIC PROBLEM - see
  % Leveque 1992, Birkhauser, Ch 6 for details!

  % Requires pulling rhol, ul, pl, but this can be done using same index
  % (ii) as the avg. states i.e. p(ii) == left pressure state for roe
  % averaged interafce (ii). P,r,u are at nodes, roe avg is at interface
  rl=r(ii); pl=p(ii); ul=u(ii); htl=ht(ii);

  LFLX = [rl*ul rl*ul^2+pl rl*htl*ul]';
%   RFLX = [rr*ur rr*ur^2+pr rr*htr*ur]';

  for ll=1:3
    lwv=lams(ll);
    lambda=min(0,lwv);
    rvec=RMAT(:,ll);
    mag=WVMAG(ll);
    % Flux is linear combination of eigenvectors and flux from left state
    % This version of scheme propagates information from left to right
    % across domain (i.e. in positive i-direction)
    LFLX=LFLX+lambda*mag*rvec;
  end

  FLXIP = LFLX;

  % Left flux at i^th+1 cell is Right flux of i^th cell
  FIP(:,ii) = FLXIP;
  FIM(:,ii+1) = FLXIP;
end

% Dirichelt BC hard fix for flux function (no net flux at domain endpoints)
FIM(:,1)=FIP(:,1);
FIP(:,nx)=FIM(:,nx);

% Transpose b/c Matlab cant understand 1D arrays well...
FIP=FIP';
FIM=FIM';

end