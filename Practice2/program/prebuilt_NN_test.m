%% TEST FRAMEWORK
clear; close all; clc;
% load predefined NN
data = load('5nets.mat');
nets = {data.net1, data.net2, data.net3, data.net4, data.net5};
% sizes of IO vectors
x_len = 31;
u_len = length(data.x1(1,:) - x_len);
% test case
x_practical = data.x1(1,1:x_len);
u_practical = data.x1(1,x_len+1:end-1);
y_practical = data.x1(1,end);
[y u] = finder(nets, x_practical, u_len);
% results
fprintf('Number of bed/days: %d\n', y);
fprintf('Control vector: \n');
for i = 1:u_len
    fprintf('%d ', u(i));
end
% error
fprintf('\n\nAbsolute error on y: %d\n', abs(y - y_practical));
fprintf('Relative error on y: %d\n', abs(y - y_practical) / y_practical);
