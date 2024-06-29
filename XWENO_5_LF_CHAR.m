function [FP1,FM1] = XWENO_5_LF_CHAR(Q,g)
%%% Uses WENO5 to compute F_j+1/2 at each grid point
%-- Also need to account for F_j-1/2

nx=size(Q,1);
neqs=size(Q,2);
FP1=zeros(nx,neqs);
VPOSRECON=zeros(neqs,1);
VMINRECON=zeros(neqs,1);
VLX=zeros(neqs,1);


%%% Needed fluxes for stencils
fm2=zeros(neqs,1);
fm1=zeros(neqs,1);
fc0=zeros(neqs,1);
fp1=zeros(neqs,1);
fp2=zeros(neqs,1);
fp3=zeros(neqs,1);

% Characteristic flux variables
vm2=zeros(neqs,1);
vm1=zeros(neqs,1);
vc0=zeros(neqs,1);
vp1=zeros(neqs,1);
vp2=zeros(neqs,1);
vp3=zeros(neqs,1);

qm2=zeros(neqs,1);
qm1=zeros(neqs,1);
qc0=zeros(neqs,1);
qp1=zeros(neqs,1);
qp2=zeros(neqs,1);
qp3=zeros(neqs,1);

%%% Get flux function from state variables
FLX=Q2FLUX(Q,g);

awvp=zeros(nx,1);
ewv=zeros(nx,1);
awvm=zeros(nx,1);

%%% Get eigenvalues for flux splitting
for ii=1:nx
  myq=Q(ii,:);
  [wv1,wv2,wv3]=EULERWAVES(myq,g);
  awvp(ii)=wv1;
  ewv(ii)=wv2;
  awvm(ii)=wv3;
end

%%% Global flux splitting values
alpha1=max(abs(awvp));
alpha2=max(abs(ewv)); % <--- order needs to be consistent with order of eigenvectors
alpha3=max(abs(awvm));


for ii=1:nx

  if ii<3 % <--- Reduce to no flux for ease
    fm2(:)=0;
    fm1(:)=0;
    fc0(:)=0;
    fp1(:)=0;
    fp2(:)=0;
    fp3(:)=0;

    qm2(:)=0;
    qm1(:)=0;
    qc0(:)=0;
    qp1(:)=0;
    qp2(:)=0;
    qp3(:)=0;


  elseif ii>nx-3 % <--- Reduce to no flux for ease
    fm2(:)=0;
    fm1(:)=0;
    fc0(:)=0;
    fp1(:)=0;
    fp2(:)=0;
    fp3(:)=0;

    qm2(:)=0;
    qm1(:)=0;
    qc0(:)=0;
    qp1(:)=0;
    qp2(:)=0;
    qp3(:)=0;

  else
    fm2(:)=FLX(ii-2,:);
    fm1(:)=FLX(ii-1,:);
    fc0(:)=FLX(ii,:);
    fp1(:)=FLX(ii+1,:);
    fp2(:)=FLX(ii+2,:);
    fp3(:)=FLX(ii+3,:);

    qm2=Q(ii-2,:)';
    qm1=Q(ii-1,:)';
    qc0=Q(ii,:)';
    qp1=Q(ii+1,:)';
    qp2=Q(ii+2,:)';
    qp3=Q(ii+3,:)';

  end

  %%% Once all candidate fluxes for stencils are obtained, do
  %%% characteristex variable mapping

  qc0=Q(ii,:)';
  if ii==nx
    qp1=qc0;
  else
    qp1=Q(ii+1,:)';
  end

  [RMAT,LMAT,LAMS] = PRIM2CHAR(qc0,qp1,g,ii);

  vm2=LMAT*fm2;
  vm1=LMAT*fm1;
  vc0=LMAT*fc0;
  vp1=LMAT*fp1;
  vp2=LMAT*fp2;
  vp3=LMAT*fp3; % <--- order needs to be consistent with order of eigenvectors

  VMAT=[vm2 vm1 vc0 vp1 vp2 vp3];

  wm2=LMAT*qm2;
  wm1=LMAT*qm1;
  wc0=LMAT*qc0;
  wp1=LMAT*qp1;
  wp2=LMAT*qp2;
  wp3=LMAT*qp3;

  WMAT=[wm2 wm1 wc0 wp1 wp2 wp3];

  [VPLS,VMNS] = CHARSPLIT(VMAT,WMAT,alpha1,alpha2,alpha3);


  for ll=1:neqs

    POSIN=VPLS(ll,:);

    VPOSRECON(ll) = WENO5_RIGHTWAVES(POSIN,ii);


    NEGIN=VMNS(ll,:);


    VMINRECON(ll) = WENO5_LEFTWAVES(NEGIN,ii);

  end

  % Combine
  VLX = VPOSRECON + VMINRECON;

  % Convert back to physical flux
  PHYSFLX=RMAT*VLX;
  FP1(ii,:)=PHYSFLX;

  if ii==101
    db=1;
  end

end



FM1 = circshift(FP1,1,1);
% FM1=FP1;

FM1(1:4,:) = FP1(1:4,:);
FM1(nx-2:nx,:) = FP1(nx-2:nx,:);

RES=FP1-FM1;
db=1;






end