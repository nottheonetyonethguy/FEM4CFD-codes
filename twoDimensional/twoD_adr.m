%% 2D steady state advection diffusion

close all
clear

x = linspace(0, 1, 10);
y = linspace(0, 1, 10);

[X, Y] = meshgrid(x, y);

Z = zeros(10, 10);
% surf(X, Y, Z)

%% testing fops
fid = fopen("test.txt", "r")
line = fgets(fid)
linestr = strsplit(line, ",")
