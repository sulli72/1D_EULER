function flxout = WENO5_LEFTWAVES(flxvecin,ii)

fm2 = flxvecin(1);
fm1 = flxvecin(2);
fc0 = flxvecin(3);
fp1 = flxvecin(4);
fp2 = flxvecin(5);
fp3 = flxvecin(6);

% Linear weights for O(dx^3) stencil
d0=3/10; %3/10;
d1=3/5;
d2=1/10; %1/10;
eps=1E-8; % div. by 0 fix


% Weights for ENO Stencils
c0j = [1/3 5/6 -1/6]; % I+1,I,I-1 stencil (r=0 stencil)
c1j = [-1/6 5/6 1/3]; % I+2,I+1,I stencil (r=1 stencil)
c2j = [1/3 -7/6 11/6]; % I+3,I+2,I+1 stencil (r=2 stencil)

S0 = [fp1 fc0 fm1];
S1 = [fp2 fp1 fc0];
S2 = [fp3 fp2 fp1];

%   fvup = [fm1 fc0 fp1];
%   fvcn = [fc0 fp1 fp2];
%   fvdn = [fp1 fp2 fp3];

fr0 = c0j*S0';
fr1 = c1j*S1';
fr2 = c2j*S2';

%   fintm = cr0*fvup';
%   fintc = cr1*fvcn';
%   fintp = cr2*fvdn';

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

% if ii==101;
% fprintf(1,'LEFT, w0 = %f, w1 = %f, w2 = %f\n',w0,w1,w2);
% end
% if ii==101;
% fprintf(1,'LEFT, B0 = %f, B1 = %f, B2 = %f\n',B0,B1,B2);
% end

ws=w0+w1+w2;
if abs(1-ws) > 0.00001
  fprintf(1,'Bad Weight Combo, Wsum = %f\n',ws)
end

flxout = w0*fr0 + w1*fr1 + w2*fr2;


end