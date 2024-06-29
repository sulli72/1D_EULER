function CHARPLOT(VAR,xv,tv,varname)
% Subroutine to plot the primitive variables vs time 
% --- helps visualize the fluid waves in a given problem

PLOTMAT=VAR';
tstr=strcat('Characteristic Field for~',varname);
[X,T] = meshgrid(xv,tv);
figure();
% p=pcolor(X,T,PLOTMAT);
[p,h]=contourf(X,T,PLOTMAT,64);
set(h,'EdgeColor','none');
colormap jet
colorbar 
title(tstr);

end