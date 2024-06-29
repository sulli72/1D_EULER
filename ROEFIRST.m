function QNP1 = ROEFIRST(Q,dx,dt,g)
% 1st Order Euler Explicit Time Marching

% Inputs:
% Q == solution at current timestep
% dx == grid spacing
% dt == time step
% g == gamma

% Outputs:
% QNP1 == solution after RK3 algorithm marches forward in time

%%% General Conservative formulation: 
%  UNP1(j) = UN(j) - dt/dx ( F(i+1/2) - F(i-1/2) ) 

% Only use first order Euler scheme for time integration.
% Explicit scheme used for ease of implementation.
% Not much to be gained by using higher order time marching schemes here,
% since Roe scheme is a highly dissipative 1st order upwind spatial scheme

[FIP,FIM] = ROEFLUX(Q,g); 
QNP1 = Q - dt/dx*(FIP-FIM);
end