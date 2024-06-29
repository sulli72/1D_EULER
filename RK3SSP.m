function QNP1 = RK3SSP(Q,dx,dt,g,fluxopt)
% 3rd Order TVD Runge-Kutta Explicit Time Marching
% Follows scheme developed and proved by Gottlieb and Shu:
% "TOTAL VARIATION DIMINISHING RUNGE-KUTTA SCHEMES - 1998 
% (Journal: Mathematics of Computation)

% Inputs:
% Q == solution at current timestep
% dx == grid spacing
% dt == time step
% g == gamma
% fluxopt == selected flux scheme

% Outputs:
% QNP1 == solution after RK3 algorithm marches forward in time

%%% General Conservative formulation:
%  UNP1(j) = UN(j) - dt/dx ( F(i+1/2) - F(i-1/2) )

if fluxopt==2

  % Weno 5 w/ Lax-Friedrichs, component wise flux reconstruction
  [FIP,FIM] = XWENOFLUX(Q,g);
  Q1 = Q - dt/dx*(FIP - FIM);

  [FIP,FIM] = XWENOFLUX(Q1,g);
  Q2 = 3/4*Q + 1/4*(Q1 - dt/dx*(FIP-FIM));

  [FIP,FIM] = XWENOFLUX(Q2,g);
  QNP1 = 1/3*(Q) + 2/3*(Q2 - dt/dx*(FIP - FIM));

elseif fluxopt==3

  % Weno 5 w/ Lax-Friedrichs, characteristic wise flux reconstruction
  [FIP,FIM] = XWENO_5_LF_CHAR(Q,g);
  Q1 = Q - dt/dx*(FIP - FIM);

  [FIP,FIM] = XWENO_5_LF_CHAR(Q1,g);
  Q2 = 3/4*Q + 1/4*(Q1 - dt/dx*(FIP-FIM));

  [FIP,FIM] = XWENO_5_LF_CHAR(Q2,g);
  QNP1 = 1/3*(Q) + 2/3*(Q2 - dt/dx*(FIP - FIM));

elseif fluxopt==4

  % Weno 5 w/ Roe Formulation, characteristic wise flux reconstruction
  [FIP,FIM] = XWENORF(Q,g);
  Q1 = Q - dt/dx*(FIP - FIM);

  [FIP,FIM] = XWENORF(Q1,g);
  Q2 = 3/4*Q + 1/4*(Q1 - dt/dx*(FIP-FIM));

  [FIP,FIM] = XWENORF(Q2,g);
  QNP1 = 1/3*(Q) + 2/3*(Q2 - dt/dx*(FIP - FIM));

elseif fluxopt==5

  % Weno 6 w/ Lax-Friedrichs, characteristic wise flux reconstruction
  [FIP,FIM] = XWENO_6_LF_CHAR(Q,g);
  Q1 = Q - dt/dx*(FIP - FIM);

  [FIP,FIM] = XWENO_6_LF_CHAR(Q1,g);
  Q2 = 3/4*Q + 1/4*(Q1 - dt/dx*(FIP-FIM));

  [FIP,FIM] = XWENO_6_LF_CHAR(Q2,g);
  QNP1 = 1/3*(Q) + 2/3*(Q2 - dt/dx*(FIP - FIM));

elseif fluxopt==6

  % Weno 6 w/ Roe formulation, characteristic reconstruction
  [FIP,FIM] = CU6WENO(Q,g);
  Q1 = Q - dt/dx*(FIP - FIM);

  [FIP,FIM] = CU6WENO(Q1,g);
  Q2 = 3/4*Q + 1/4*(Q1 - dt/dx*(FIP-FIM));

  [FIP,FIM] = CU6WENO(Q2,g);
  QNP1 = 1/3*(Q) + 2/3*(Q2 - dt/dx*(FIP - FIM));

elseif fluxopt>6
  fprintf(1,'Invalid flux scheme choice\n')
end

end