clear;
clc;
close all;

data = importdata("Laplace_sol.txt");
x = data(:,1)';
y = data(:,2)';
z = data(:,3)';

figure;
plot3(x,y,z,'o');
xlabel('x');
ylabel('y');
zlabel('z');

