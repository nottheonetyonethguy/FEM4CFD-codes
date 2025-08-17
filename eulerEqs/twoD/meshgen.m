% function [] = meshgen(ndim, ndof, xL, xR, yL, yR, n_elem, matdata, k, ux, uy, elemdata, e1, e2, npelem, nforce)
clc; clear; % close all;
ndim = 2; % no of dimensions
ndof = 4; % no of dofs
n_elem_x = 20; % no of elements x direction
n_elem_y = 4; % no of elements y direction
nnode = (n_elem_x + 1) * (n_elem_y + 1);

% domain boundaries
xL = 0; xR = 1;
yL = 0; yR = 0.1;

% matdata
matdata = 1; gamma = 1.4; f = 0;
rho_l = 1.0; rho_r = 0.125;
u_l = 0; u_r = 0;
v_l = 0; v_r = 0;
p_l = 1; p_r = 0.1;

% elemdata
elemdata = 1;
e1 = 0.0;
e2 = 0.0;

npelem = 4;
nelem = n_elem_x * n_elem_y;

% boundary conditions
nforce = 0;
% nbound = 100;
% number of nodes in boundary = n_elem * 4 (sides)
% nbound = (n_elem_x + 1) * (n_elem_y + 1);
nbound = (n_elem_y + 1) * 2; % left most and right most sides

%% open new file to write data
% delete(data.inp)
fileID = fopen('data0.inp', 'a+');
fprintf(fileID, 'ndim,%d\n', ndim);
fprintf(fileID, 'ndof,%d\n', ndof);
fprintf(fileID, 'nnode,%d\n', nnode);

%% node printing

x = linspace(xL, xR, n_elem_x + 1);
y = linspace(yL, yR, n_elem_y + 1);

[X, Y] = meshgrid(x, y);
Z = zeros(n_elem_x + 1, n_elem_y + 1);

surf(X, Y, Z');

node_coords = zeros(nnode, 2);

for i = 0:n_elem_y
	x_coords = X(i + 1, :)';
	y_coords = Y(i + 1, :)';
	ki = 1;
	for j = (((n_elem_x + 1) * i) + 1):((n_elem_x + 1) * (i + 1))
		node_coords(j, 1) = x_coords(ki);
		node_coords(j, 2) = y_coords(i + 1);
		ki = ki + 1;
		fprintf(fileID, '%d,%f,%f\n', j, node_coords(j, 1), node_coords(j, 2));
	end
end

%% matdata

fprintf(fileID, 'matdata,%d\n', matdata);

for i = 1:matdata
	fprintf(fileID, '%d,%.8f,%f\n', i, gamma, f);
end

%% elemdata

fprintf(fileID, 'elemdata,%d\n', elemdata);

for i = 1:elemdata
	fprintf(fileID, '%d,%f,%f\n', i, e1, e2);
end

fprintf(fileID, 'npelem,%d\n', npelem);
fprintf(fileID, 'nelem,%d\n', nelem);

% elem connectivity

% generate conn_table
conn_table = zeros(n_elem_y + 1, n_elem_x + 1);
xx = n_elem_x + 1;

for i = 0:n_elem_y
	
	for j = 1:(n_elem_x + 1)
		conn_table(i + 1, j) = xx * i + j;
	end
	
end

% print element connectivity table
xx = 1;

for i = 0:n_elem_y - 1
	l1 = conn_table(i + 1, :)';
	l2 = conn_table(i + 2, :)';
	
	for j = 0:n_elem_x - 1
		fprintf(fileID, '%d,%d,%d,%d,%d,%d,%d\n', xx, matdata, elemdata, l1(j + 1), l1(j + 2), l2(j + 2), l2(j + 1));
		xx = xx + 1;
	end
	
end

%% setting initial values
fprintf(fileID, "ninit,%d\n", nnode);
x_coords = X(1, :); % x coordinates to check the initial condition
for i = 1:(n_elem_x + 1)
	node_num = conn_table(:, i);
	if(x_coords(i) <= 0.5)
		rho = rho_l;
		u = u_l; v = v_l;
		p = p_l;
	else
		rho = rho_r;
		u = u_r; v = v_r;
		p = p_r;
	end
	for j = 1:max(size(node_num))
		fprintf(fileID, '%d,%d,%f,%f,%f,%f\n', node_num(j), elemdata, rho, u, v, p);
	end
end

%% boundary layers
fprintf(fileID, "nbound,%d\n", nbound);
left = conn_table(:, 1); right = conn_table(:, end);
node_nums = [left right]; node_nums = node_nums(:);
for i = 1:nbound
	if (i <= 5)
		rho = rho_l; u = u_l; v = v_l; p = p_l;
	else
		rho = rho_r; u = u_r; v = v_r; p = p_r;
	end
	fprintf(fileID, '%d,%d,%f,%f,%f,%f\n', node_nums(i), elemdata, rho, u, v, p);
end

%% forcing term
fprintf(fileID, "nforce,%d\n", nforce);

% end
