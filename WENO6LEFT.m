function flxout = WENO6LEFT(flxvecin,ii)

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
c0j = [1/3 5/6 -1/6]; % I+1,I,I-1 stencil (r=0 stencil)
c1j = [-1/6 5/6 1/3]; % I+2,I+1,I stencil (r=1 stencil)
c2j = [1/3 -7/6 11/6]; % I+3,I+2,I+1 stencil (r=2 stencil)
cdwj = [11/6 -7/6 1/3]; % I,I-1,I-2 stencil (added downwind stencil)


S0 = [fp1 fc0 fm1];
S1 = [fp2 fp1 fc0];
S2 = [fp3 fp2 fp1];
SDW = [fc0 fm1 fm2];

fr0 = c0j*S0';
fr1 = c1j*S1';
fr2 = c2j*S2';
frdw = cdwj*SDW';

% Compute smoothness indicator polynomials
B0=13/12*(fp1 - 2*fc0 + fm1)^2 + 1/4*(3*fp1 - 4*fc0 + fm1)^2; %I-1,I,I+1 stencil (measures 1st, 2nd derivs using FD of stencil cell avgs)
B1=13/12*(fp2 - 2*fp1 + fc0)^2 + 1/4*(fp2 - fc0)^2; %I,I+1,I+2 stencil
B2=13/12*(fp3 - 2*fp2 + fp1)^2 + 1/4*(fp3 -4*fp2 + 3*fp1)^2; %I+1,I+2,I+3 stencil

%%% ==== WEIGHTS FROM BREHM PAPER === %%%
B6 = 1/120960 *( 271779*fp3^2 + fp3*(-2380800*fp2 + 4086352*fp1 - 3462252*fc0 + 1458762*fm1 - 245620*fm2) ...
  + fp2*(5653317*fp2 - 20427884*fp1 + 17905032*fc0 - 7727988*fm1 + 1325006*fm2) ...
  + fp1*(19510972*fp1 - 35817664*fc0 + 15929912*fm1 - 2792660*fm2) ...
  + fc0*(17195652*fc0 - 15880404*fm1 + 2863984*fm2) ...
  + fm1*(3824847*fm1 - 1429976*fm2) + 139633*fm2^2);



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

% if ii==101;
% fprintf(1,'LEFT, B0 = %f, B1 = %f, B2 = %f, B6 = %f\n',B0,B1,B2,B6);
% end

% if ii==101;
% fprintf(1,'%i ~~~ %f ~~~ %f  ~~~ %f  ~~~ %f\n',CFC,Tau6/(B0+eps),Tau6/(B1+eps),Tau6/(B2+eps),Tau6/(B3+eps))
% 
% db=1;
% 
% end

%   for jj=1:4
%     myB=wvec(jj);
%     if myw<-0.01
%       fprintf(1,'Negative left weight\nii = %i, j = %i wt = %f\n\n',ii,jj,myw)
%     end
%   end

ws=w0+w1+w2+wdw;
if abs(1-ws) > 0.01
  fprintf(1,'Bad Weight Combo, Wsum = %f\n',ws)
end

flxout =  w0*fr0 + w1*fr1 + w2*fr2 + wdw*frdw;

end