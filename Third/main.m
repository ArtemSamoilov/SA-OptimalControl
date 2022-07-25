% Samoilov Artem, 12.09.2021
% Task 1 for System Analysis Practicum on Optimal Control
% Unauthorizd use of creator's intellectual property is strictly forbidden
%%
clc
clear
alpha = 20;%2,3,4,5
t = 0.5;
[X,Y,x_l,y_l] = reachset(alpha, t);
hold on
%axis([-2 1.5 -2.2 2]);
plot(X,Y,'Color','b','LineWidth',3);
plot(x_l,y_l,'Color','r','LineWidth',3);
xlabel('x_1');
ylabel('x_2');
hold off

t_1 = 0.2;
t_2 = 0.5;
filename = 'name.avi';
N = 10;
%reachsetdyn(alpha,t_1,t_2,N,filename);





