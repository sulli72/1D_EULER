function PRIMPLOT(r,u,p,xv,g)
% Subroutine to plot the primitive variables
% NOTE: 
% u = velocity
% p = pressure
% r = density

figure(1)
plot(xv,r,'c-',xv,u,'r-',xv,p,'b-','LineWidth',2);
legend('$\rho$','$u$','$P$','Location','best')
grid on
xlabel('$x$')
ylabel('$U^*$')


end