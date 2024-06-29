%% Main driver code for solving 1D Euler system
% The compressible gas-dynamics equation are cast in conservative form
% q_t + F(q)_x = 0
%
% q = [r ru E];
% q is the vector of conserved variables
% r == density
% u == velocity
% E == total energy ( E = re, e = p/(r(gamma - 1)) + 0.5*u^2 );
% P == pressure 

% Multiple spatial discretizations for hyperbolic system of eqns the can be selected
%
% Option 1: 1st Order Roe scheme -- the classical approximate Riemann
% solver of Prof. Philip Roe
%
% Option 2: Component-wise 5th order WENO scheme of Jang and Shu, with Lax-Friedrichs flux
% splitting for stability
%
% Option 3: Characteristic-wise 5th order WENO scheme with Lax-Friedrichs
% flux splitting
%
% Option 4: Characteristic-wise 5th order WENO scheme of Jang and Shu, 'Roe
% Formulation'
%
% Option 5: 6th order centered WENO scheme of Hu and Adams
% (with Lax-Friedrichs flux splitting for stability and Roe averaging for the characteristic mapping)
%
% Option 6: 6th order centered WENO scheme of Hu and Adams
% (with Weno-Roe formulation for the upwinding and characteristic decomposition)
% -- this approach requires special modification of non-linear WENO weights
% to improve stability
%
% SOME USEFUL SOURCES:
% Roe Scheme:
%  "Approximate Riemann solvers, parameter vectors, and difference schemes"
%  by P. Roe, 1981 (original paper, Journal of Computational Physics)
%  "Copmutational Gasdynamics" by C. Laney (textbook)
%
% WENO5 Scheme
% "Essentially Non-Oscillatory and Weighted Essentially Non-Oscillatory
% Schemes for Hyperbolic Conservation Laws" by C.W. Shu (ICASE/NASA 1997 Report)
%
% "A comparison of higher-order finite-difference shock capturing schemes"
%  by Brehm et al (Computers and Fluids 2015)
%
% WENO6 Scheme
% "A comparison of higher-order finite-difference shock capturing schemes"
%  by Brehm et al (Computers and Fluids 2015)
% 
% "An adaptive central-upwind weighted essentially non-oscillatory scheme"
% by Hu et al (Journal of Computational Physics 2010)

clearvars; clc; close all;
set(0,'defaultFigureRenderer','painters')
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
%%% ===== End of Header ===== %%%

%% Mesh and Simulation Parameters
% Number of grid points
nx=200;

% Ratio of specific heats for gas of choice - g == 1.4 for air
g=1.4;

% Choice of flux scheme (see header for details)
flxopt=3; % <-- should be in the range 1-6

% Choose whether to visualize characteristic waves (0=no,1=yes)
iplotchar=1;


%% CASE OPTIONS

% Uncomment the case you want to use
% CASE='SOD'; 
% CASE='SHU'; 
% CASE='STEADY'; 
% CASE='LAX'; 
CASE='123'; 

% Set to 1 to visualize solution as it marches in time
% Set to 0 to NOT visualze during solution, but only plot solution after
% final timestep
imvie=0;

% Alloate initial condition
[Q,xv,tf]=ICMAKER(CASE,nx,g);
dx=xv(2)-xv(1);
[r,u,p] = Q2PRIM(Q,g);

% Plots initial condition
PRIMPLOT(r,u,p,xv,g)

% Speed of sound, used for flux splitting and CFL constraint
c=(g*p./r).^0.5;

% fprintf(1,'PLOTTING INITIAL CONDITIONS FOR CHOSE PROBLEM\n')
% cont=input('Enter any number to continue\n');

%% Time Marching Parameters
CFL=.6;              % CFL number
wvsp=max(abs(u)+c);  % Maximum characteristic wavespeed in initial condition
dt=CFL*dx/wvsp;      % Timestep size -- selected for stability
nsteps=floor(tf/dt); % Determine number of timesteps


tic %<-- meaure run time
t=0;
it=0;
maxits=2E3; %<-- set upper bound on iterations if needed

while t<tf && it<maxits

  % Call subroutine used for finding maximum wavespeed in solution on the fly
  wvsp=WAVESPEED(Q,g);
  % Update dt on the fly to maintain stability
  dt=CFL*dx/wvsp;

  %%% --- 1st Order Roe Scheme w/ 1st order time integration --- %%%
  if flxopt==1
    Q = ROEFIRST(Q,dx,dt,g);   % <--- Spatial scheme option 1
  else
  %%% --- 3rd Order time integration for higher order spatial schemes --- %%%
    Q = RK3SSP(Q,dx,dt,g,flxopt);
  end


  % Upate command window
  fprintf(1,'Iteration %i\n',it)
  it=it+1;
  t=t+dt;

  
%% These will r(t), p(t), u(t) fields will be used to visualize characteristic waves if desired
  if iplotchar==1
    entry=it;
    tv(it)=t;
    [r,u,p]=Q2PRIM(Q,g);
    RMAT(:,entry)=r;
    UMAT(:,entry)=u;
    PMAT(:,entry)=p;
  end

  if imvie==1
    [r,u,p] = Q2PRIM(Q,g);
    PRIMPLOT(r,u,p,xv,g)
  end

end
toc % <-- gives user info about simulation run time


%% Plotting Solution
[r,u,p]=Q2PRIM(Q,g);
PRIMPLOT(r,u,p,xv,g)

%%% If SOD or SHU cases are selected, the 'exact' solution
%%% will be plotted once the numerical solution is completed
%%% Note: the exact solution is computed with the WENO5 characteristic
%%% scheme using 4000 grid points

switch CASE
  case 'SHU'
    exsoln=load('SHU_RHO_PROFILE_EXACT.mat');
    xex=exsoln.xv;
    rex=exsoln.r;
    RHOPLOT(r,xv,rex,xex);
end

switch CASE
  case 'SOD'
    exsoln=load('SOD_PROFILES_EXACT.mat');
    xex=exsoln.xv;
    rex=exsoln.r;
    pex=exsoln.p;
    uex=exsoln.u;
    RHOPLOT(r,xv,rex,xex);
end


%% Plotting 'Characteristics'

%%% PLOT PRIMITIVE VARIABLE DISTRIBUTIONS TO SEE CHARACTERISTIC WAVE
%%% FRONTS
if iplotchar==1
  CHARPLOT(RMAT,xv,tv,'Denstiy');
  CHARPLOT(UMAT,xv,tv,'Velocity');
  CHARPLOT(PMAT,xv,tv,'Pressure');
end
%% Mach number
CMAT = (1.4*PMAT./RMAT).^0.5;
MACHMAT = UMAT./CMAT;
MACHMAX = max(max(MACHMAT));
MACHMEAN = mean(MACHMAT,1);

%%% End of driver script