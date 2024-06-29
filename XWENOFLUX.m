function [FP1,FM1] = XWENOFLUX(Q,g)
%%% Uses WENO5 to compute F_j+1/2 at each grid point
%-- Also need to account for F_j-1/2

nx=size(Q,1);
neqs=size(Q,2);
FP1PLS=zeros(nx,neqs);
FP1MNS=zeros(nx,neqs);

% Get flux functions from solution i.e, flux = f(q)
FLX=Q2FLUX(Q,g);


%%% Split the flux functions f(q) into f(q) = (f(q)+)  +  (f(q)-)
[PLSFLXES,MNSFLXES] = EULERFLUXSPLIT(FLX,Q,g);
%%% From here, get do WENO on each of the 3 f+, f- vectors....

%%% ---- WENO5 for f^+ functions ---- %%%

% LINEAR weights for O(dx^3) stencil (3rd order ENO stencil)
d0=3/10; %3/10;
d1=3/5;
d2=1/10; %1/10;

eps=1E-8; % used for div. by 0 fix


%%% WENO ON FPLUS IN COMPONENT FORM
for ff=1:neqs
  fpls=PLSFLXES(:,ff);

  for ii=1:nx
    % OBTAIN FLUXES AT NEEDED POINTS FOR GIVEN I POINT
    if ii<=3 % <--- Reduce to no flux for ease
      % Do nothing at boundaries
      fm2=0;
      fm1=0;
      fc0=0;
      fp1=0;
      fp2=0;

    elseif ii>=nx-2 % <--- Reduce to no flux for ease
      % Do nothing at boundaries
      fm2=0;
      fm1=0;
      fc0=0;
      fp1=0;
      fp2=0;

    else
      % Get needed fluxes
      fm2=fpls(ii-2);
      fm1=fpls(ii-1);
      fc0=fpls(ii);
      fp1=fpls(ii+1);
      fp2=fpls(ii+2);

    end

    % Weights for ENO Stencils 
    % This is why this class of schemes are called 'Weighted ENO'
    c0j = [1/3 5/6 -1/6]; % I,I+1,I+2 stencil (r=0 stencil)
    c1j = [-1/6 5/6 1/3]; % I-1,I,I+1 stencil (r=1 stencil)
    c2j = [1/3 -7/6 11/6]; % I-2,I-1,I stencil (r=2 stencil)

    % Build stencils as vectors
    S0 = [fc0 fp1 fp2];
    S1 = [fm1 fc0 fp1];
    S2 = [fm2 fm1 fc0];

    % Put ENO weights onto stencil
    fr0 = c0j*S0';
    fr1 = c1j*S1';
    fr2 = c2j*S2';

    % Compute smoothness indicator polynomials
    B0=13/12*(fc0 - 2*fp1 + fp2)^2 + 1/4*(3*fc0 - 4*fp1 + fp2)^2; %I,I-1,I-2 stencil (measures 1st, 2nd derivs using FD of stencil cell avgs)
    B1=13/12*(fm1 - 2*fc0 + fp1)^2 + 1/4*(fm1 - fp1)^2; %I,I-1,I+1 stencil
    B2=13/12*(fm2 - 2*fm1 + fc0)^2 + 1/4*(fm2 - 4*fm1 + 3*fc0)^2; %I,I+1,I+2 stencil

    % Compute alpha values (needed for non-linear weights)
    a0 = d0/(eps+B0)^2;
    a1 = d1/(eps+B1)^2;
    a2 = d2/(eps+B2)^2;
    as = a0 + a1 + a2; %sum of alphas so non-linear weights normalized to add up to one

    % Compute Non-Linear Weights
    w0 = a0/as;
    w1 = a1/as;
    w2 = a2/as;

    ws=w0+w1+w2;
    % Check to see if weights sum to 1
    if abs(1-ws) > 0.00001
      fprintf(1,'Bad Weight Combo, Wsum = %f\n',ws)
    end

    % Get weighted reconstruction of f+ flux at given point
    FP1PLS(ii,ff) =  w0*fr0 + w1*fr1 + w2*fr2;
  end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ------- FLUX FOR LEFT RUNNING WAVES ----- %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% WENO using appropriate stencils for the flux from the left running wave, f-
%%% Note: the stencils for f- are the 'reverse' of the stencils for the
%%% f+ waves, to ensure upwinding for the left running waves

for ff=1:neqs
  fmns=MNSFLXES(:,ff);

  % Loop over x-direction
  for ii=1:nx 

    % OBTAIN FLUXES AT NEEDED POINTS FOR GIVEN I POINT
    if ii<=3 
      % Do nothing at boundaries
      fm1=0;
      fc0=0;
      fp1=0;
      fp2=0;
      fp3=0;

    elseif ii>=nx-2
      % Do nothing at boundaries
      fm1=0;
      fc0=0;
      fp1=0;
      fp2=0;
      fp3=0;

    else
      % Get flux at pts [ i-1 ... i+3]
      fm1=fmns(ii-1);
      fc0=fmns(ii);
      fp1=fmns(ii+1);
      fp2=fmns(ii+2);
      fp3=fmns(ii+3);
    end

    % Weights for ENO Stencils
    c0j = [1/3 5/6 -1/6]; % I+1,I,I-1 stencil (r=0 stencil)
    c1j = [-1/6 5/6 1/3]; % I+2,I+1,I stencil (r=1 stencil)
    c2j = [1/3 -7/6 11/6]; % I+3,I+2,I+1 stencil (r=2 stencil)

    % Build stencils as vectors
    S0 = [fp1 fc0 fm1];
    S1 = [fp2 fp1 fc0];
    S2 = [fp3 fp2 fp1];

    % Get flux for each stencil
    fr0 = c0j*S0';
    fr1 = c1j*S1';
    fr2 = c2j*S2';

    % Compute smoothness indicator polynomials
    B0=13/12*(fp1 - 2*fc0 + fm1)^2 + 1/4*(3*fp1 - 4*fc0 + fm1)^2; %I-1,I,I+1 stencil (measures 1st, 2nd derivs using FD of stencil cell avgs)
    B1=13/12*(fp2 - 2*fp1 + fc0)^2 + 1/4*(fp2 - fc0)^2; %I,I+1,I+2 stencil
    B2=13/12*(fp3 - 2*fp2 + fp1)^2 + 1/4*(fp3 -4*fp2 + 3*fp1)^2; %I+1,I+2,I+3 stencil


    % Compute alpha values (needed for non-linear weights)
    a0 = d0/(eps+B0)^2;
    a1 = d1/(eps+B1)^2;
    a2 = d2/(eps+B2)^2;
    as = a0 + a1 + a2;

    % Compute Non-Linear Weights
    w0 = a0/as;
    w1 = a1/as;
    w2 = a2/as;

    ws=w0+w1+w2;
    if abs(1-ws) > 0.00001
      fprintf(1,'Bad Weight Combo, Wsum = %f\n',ws)
    end

    % Get weighted reconstruction of f- flux at given point
    FP1MNS(ii,ff) = w0*fr0 + w1*fr1 + w2*fr2;

  end

end

% Build flux at i+1/2 points
FP1 = FP1PLS + FP1MNS;
% Flux at i-1/2 is a circular shift of flux at i+1/2... draw it out to see.
FM1 = circshift(FP1,1,1);

% HACK -- SETTING FLUXES TO ZERO AT THE BOUNDARIES SO SOLN DOESNT UPDATE AT THESE
% POINTS... THIS IS JUST FOR EASE OF IMPLEMENTATION DO NOT DO THIS IN
% ACTUAL CODE.
FM1(1:4,:) = FP1(1:4,:);
FM1(nx-2:nx,:) = FP1(nx-2:nx,:);
end