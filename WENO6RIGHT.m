function flxout = WENO6RIGHT(flxvecin,ii)

fm2 = flxvecin(1);
fm1 = flxvecin(2);
fc0 = flxvecin(3);
fp1 = flxvecin(4);
fp2 = flxvecin(5);
fp3 = flxvecin(6);


% Linear weights for O(dx^3) stencils
d0=9/20; %3/10;
d1=9/20;
d2=1/20; %1/10;
ddw=1/20;
eps=1E-16; % div. by 0 fix
CFC=20;
pwr=1; %<--- adjusting this to be >1 minimizes contact discontinuity wiggles (WENO6-RF FORMULATION)

% Weights for ENO Stencils
c0j = [1/3 5/6 -1/6]; % I,I+1,I+2 stencil (r=0 stencil)
c1j = [-1/6 5/6 1/3]; % I-1,I,I+1 stencil (r=1 stencil)
c2j = [1/3 -7/6 11/6]; % I-2,I-1,I stencil (r=2 stencil)
cdwj = [11/6 -7/6 1/3]; % I+1,I+2,I+3 stencil (r=-1 stencil)

S0 = [fc0 fp1 fp2];
S1 = [fm1 fc0 fp1];
S2 = [fm2 fm1 fc0];
SDW = [fp1 fp2 fp3];


fr0 = c0j*S0';
fr1 = c1j*S1';
fr2 = c2j*S2';
frdw = cdwj*SDW';

% Compute smoothness indicator polynomials
B0=13/12*(fc0 - 2*fp1 + fp2)^2 + 1/4*(3*fc0 - 4*fp1 + fp2)^2; %I,I-1,I-2 stencil (measures 1st, 2nd derivs using FD of stencil cell avgs)
B1=13/12*(fm1 - 2*fc0 + fp1)^2 + 1/4*(fm1 - fp1)^2; %I,I-1,I+1 stencil
B2=13/12*(fm2 - 2*fm1 + fc0)^2 + 1/4*(fm2 - 4*fm1 + 3*fc0)^2; %I,I+1,I+2 stencil

  %%% ==== SMOOTHNESS INDICATOR FROM BREHM PAPER ==== %%%
B6 = 1/120960 *( 271779*fm2^2 + fm2*(-2380800*fm1 + 4086352*fc0 - 3462252*fp1 + 1458762*fp2 - 245620*fp3) ...
  + fm1*(5653317*fm1 - 20427884*fc0 + 17905032*fp1 - 7727988*fp2 + 1325006*fp3) ...
  + fc0*(19510972*fc0 - 35817664*fp1 + 15929912*fp2 - 2792660*fp3) ...
  + fp1*(17195652*fp1 - 15880404*fp2 + 2863984*fp3) ...
  + fp2*(3824847*fp2 -1429976*fp3) + 139633*fp3^2);
%%% ABOVE B6 GIVES O(h^6) CONVERGENCE FOR BURGERS NL EQN


%%% === PROPER FORWARD WEIGHTS FOR HU AND ADAMS PAPER
% B6 = 1/10080*(271779*fm2^2 + fm2*(2380800*fm1 + 4086352*fc0 - 3462252*fp1 + 1458762*fp2 - 245620*fp3) ...
%  + fm1*(5653317*fm1 - 20427884*fc0 + 17905032*fp1 - 7727988*fp2 + 1325006*fp3)...
%  + fc0*(19510972*fc0 - 35817664*fp1 + 15929912*fp2 -2792660*fp3) ... 
%  + fp1*(17195652*fp1 - 15880404*fp2 + 2863984*fp3) ...
%  + fp2*(3824847*fp2 - 1429976*fp3) + 139633*fp3^2); 


Tau6 = B6 - 1/6*(B0 + B2 + 4*B1);
B3 = B6;


a0 = d0*(CFC + Tau6/(B0+eps))^pwr;
a1 = d1*(CFC + Tau6/(B1+eps))^pwr;
a2 = d2*(CFC + Tau6/(B2+eps))^pwr;
adw = ddw*(CFC + Tau6/(B3+eps))^pwr;


as = a0 + a1 + a2 + adw;

% Compute Non-Linear Weights
w0 = a0/as;
w1 = a1/as;
w2 = a2/as;
wdw = adw/as;

wvec=[w0 w1 w2 wdw];



ws=w0+w1+w2+wdw;
if abs(1-ws) > 0.01
  fprintf(1,'Bad Weight Combo, Wsum = %f\n',ws)
end

flxout =  w0*fr0 + w1*fr1 + w2*fr2 + wdw*frdw;
end

