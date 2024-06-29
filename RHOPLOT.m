function RHOPLOT(r,xv,rex,xex)
% Subroutine to compare exact and numerically computed density profiles
% used for Sod and Shu problems -- follows standard practice from CFD
% literature on scheme assessment.


nx=length(r);
tsr=strcat('$N_x =~$',num2str(nx));
figure()
plot(xex,rex,'b-',xv,r,'r--','LineWidth',2);
grid on
xlabel('$x$')
ylabel('$\rho$')
title(tsr);
legend('$\rho_{ex}$','$\rho_{num}$')

end