function [R,RU,E] = Q2CONS(Q,g)
% This function turns solution matrix Q into its conservative variables
% Q = [R, RU, E] ---> cons1 = R, cons2=RU, cons3=E

R=Q(:,1);
RU=Q(:,2);
E=Q(:,3);

end