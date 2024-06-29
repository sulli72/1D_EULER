function [FP1,FM1] = XWENORF(Q,g)
%%% Uses WENO5 to compute F_j+1/2 at each grid point
%-- Also need to account for F_j-1/2

nx=size(Q,1);
neqs=size(Q,2);

FP1=zeros(nx,neqs);

% Get flux function from state variables
FLX=Q2FLUX(Q,g);


%%% WENO ON FLUX IN CHARACTERISTIC FORM
neqs=3;
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

tmp=zeros(neqs,1);

for ii=1:nx

  if ii<3 % <--- Reduce to no flux for ease
    fm2(:)=0;
    fm1(:)=0;
    fc0(:)=0;
    fp1(:)=0;
    fp2(:)=0;
    fp3(:)=0;

  elseif ii>nx-3 % <--- Reduce to no flux for ease
    fm2(:)=0;
    fm1(:)=0;
    fc0(:)=0;
    fp1(:)=0;
    fp2(:)=0;
    fp3(:)=0;

  else
    fm2(:)=FLX(ii-2,:);
    fm1(:)=FLX(ii-1,:);
    fc0(:)=FLX(ii,:);
    fp1(:)=FLX(ii+1,:);
    fp2(:)=FLX(ii+2,:);
    fp3(:)=FLX(ii+3,:);
    db=1;
  end

  %%% Once all candidate fluxes for stencils are obtained, do
  %%% characteristex variable mapping

  qc0=Q(ii,:)';
  if ii==nx
    qp1=qc0;
  else
    qp1=Q(ii+1,:)';
  end  

  % Map to characteristic variables
  [RMAT,LMAT,LAMS] = PRIM2CHAR(qc0,qp1,g,ii);
%   LAMS=flip(LAMS);
% LAMS=circshift(LAMS,-1);

  vm2=LMAT*fm2;
  vm1=LMAT*fm1;
  vc0=LMAT*fc0;
  vp1=LMAT*fp1;
  vp2=LMAT*fp2;
  vp3=LMAT*fp3; % <--- order needs to be consistent with order of eigenvectors

  lcnt=0;
  rcnt=0;

  for ll=1:neqs
    gm2=vm2(ll);
    gm1=vm1(ll);
    gc0=vc0(ll);
    gp1=vp1(ll);
    gp2=vp2(ll);
    gp3=vp3(ll);

    gvec=[gm2 gm1 gc0 gp1 gp2 gp3];
    


    % Perform WENO on charactertistic variables
    % w/ direction of stencil set by sign of local eigenvalue (equivalent to wave propagation velocity)
    if LAMS(ll)>=0

        tmp(ll) = WENO5_RIGHTWAVES(gvec,ii);

    elseif LAMS(ll)<0

        tmp(ll) = WENO5_LEFTWAVES(gvec,ii);

    end

  end


% Convert back to physical flux
PHYSFLX=RMAT*tmp;
FP1(ii,:)=PHYSFLX;

  if ii==101
    db=1;
  end
end

FM1 = circshift(FP1,1,1);

%%% TEST PROBLEM HARD BC FIX (kills waves at boundaries)
FM1(1:4,:) = FP1(1:4,:);
FM1(nx-2:nx,:) = FP1(nx-2:nx,:);


end