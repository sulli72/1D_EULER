function [r,u,p] = Q2PRIM(Q,g)
% This function turns solution matrix Q into its primitive variables
% Q = [R, RU, E] ---> prim = [r, u, p];

R=Q(:,1);
RU=Q(:,2);
E=Q(:,3);
[r,u,p] = CONS2PRIM(R,RU,E,g);
end