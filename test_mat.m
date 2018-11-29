%tests my matrix class
clear
clc
close all

data = importdata('output.txt');
figure;
for i=1:length(data(:,1))
   plot(data(i,2:end));
   %axis([0 length(data(1,2:end)) 0 1+.4]);
   pause(.1);
end
% %plot(data(end,2:end))
