function [Q0,xv,tf] = ICMAKER(CASE,nx,g)
% Inputs: Case name (SOD, SHU, etc.), number of grid points, value of gamma
% This function builds the initial condition for a given case

% Initial Condition
nvars=3;
Q0=zeros(nx,nvars);
p=zeros(nx,1);
u=zeros(nx,1);
r=zeros(nx,1);

switch CASE
  case 'SOD'
    %%%% Sod Shock Tube IC

    x0=0;
    xf=1;
    dx=(xf-x0)/(nx-1);
    xv=x0:dx:xf;

    xdphr=interp1(xv,xv,0.5,'nearest');
    idphr=find(xv==xdphr);

    r(1:idphr)=1;
    u(1:idphr)=0;
    p(1:idphr)=1;

    r(idphr+1:nx)=.125;
    u(idphr+1:nx)=0;
    p(idphr+1:nx)=.1;

    % Convert to conservative solution form
    [R,RU,E] = PRIM2CONS(r,u,p,g);
    Q0 = [R RU E];

    tf=.2; %for Sod Problem

  case 'SHU'
    x0=-5;
    xf=5;
    dx=(xf-x0)/(nx-1);
    xv=x0:dx:xf;

    xdphr=interp1(xv,xv,-4,'nearest');
    idphr=find(xv==xdphr);

    r(1:idphr)=3.857143;
    u(1:idphr)=2.629369;
    p(1:idphr)=10.3333;

    eps=0.2;
    %     eps=0;
    r(idphr+1:nx)= 1 + eps*sin(5*xv(idphr+1:nx));
    u(idphr+1:nx)=0;
    p(idphr+1:nx)=1;

    % Convert to conservative solution form
    [R,RU,E] = PRIM2CONS(r,u,p,g);
    Q0 = [R RU E];

    tf=1.8; %for Shu-Osher problem

  case 'STEADY'

    x0=-1;
    xf=1;
    dx=(xf-x0)/(nx-1);
    xv=x0:dx:xf;

    xdphr=interp1(xv,xv,0,'nearest');
    idphr=find(xv==xdphr);

    r(1:idphr)=5.6;
    u(1:idphr)=1;
    p(1:idphr)=1;


    r(idphr+1:nx)= 14.933333;
    u(idphr+1:nx)=0.375;
    p(idphr+1:nx)=4.5;

    % Convert to conservative solution form
    [R,RU,E] = PRIM2CONS(r,u,p,g);
    Q0 = [R RU E];
    %     Q0=flip(Q0,1);

    tf=3; % arbitrary for steady Riemann problem

  case 'LAX'
    %%%% Lax Problem IC

    x0=0;
    xf=1;
    dx=(xf-x0)/(nx-1);
    xv=x0:dx:xf;

    xdphr=interp1(xv,xv,0.5,'nearest');
    idphr=find(xv==xdphr);

    r(1:idphr)=0.445;
    u(1:idphr)=0.698;
    p(1:idphr)=3.528;

    r(idphr+1:nx)=.5;
    u(idphr+1:nx)=0;
    p(idphr+1:nx)=.571;

    % Convert to conservative solution form
    [R,RU,E] = PRIM2CONS(r,u,p,g);
    Q0 = [R RU E];

    tf=.14; %for Lax Problem

  case '123'
    %%%% 123 Test Problem IC

    x0=0;
    xf=1;
    dx=(xf-x0)/(nx-1);
    xv=x0:dx:xf;

    xdphr=interp1(xv,xv,0.5,'nearest');
    idphr=find(xv==xdphr);

    r(1:idphr)=1;
    u(1:idphr)=-2;
    p(1:idphr)=.4;

    r(idphr+1:nx)=1;
    u(idphr+1:nx)=2;
    p(idphr+1:nx)=.4;

    % Convert to conservative solution form
    [R,RU,E] = PRIM2CONS(r,u,p,g);
    Q0 = [R RU E];

    tf=.1; %for Sod Problem

end

end