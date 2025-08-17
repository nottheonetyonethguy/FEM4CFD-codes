% function [] = meshgen(ndim, ndof, xL, xR, yL, yR, n_elem, matdata, k, ux, uy, elemdata, e1, e2, npelem, nforce)
clc; clear; % close all;
ndim = 2; % no of dimensions
ndof = 1; % no of dofs
n_elem = 20; % no of elements
nnode = (n_elem + 1) * (n_elem + 1);

% domain boundaries
xL = 0; xR = 1;
yL = 0; yR = 1;

% matdata
matdata = 1;
% k = 0.0001;
k = 1e-8;
% ux = 0.7071; uy = 0.7071;
% case 1, delta = 0
mod_u = 1;
delta = 0;
ux = mod_u * cosd(delta); uy = -mod_u * sind(delta);
s = 0; f = 0;

% elemdata
elemdata = 1;
e1 = 0.0;
e2 = 0.0;

npelem = 4;
nelem = n_elem ^ 2;

% boundary conditions
nforce = 0;
% nbound = 100;
% number of nodes in boundary = n_elem * 4 (sides)
nbound = n_elem * 4;

%% open new file to write data
% delete(data.inp)
fileID = fopen('data0.inp', 'a+');
fprintf(fileID, 'ndim,%d\n', ndim);
fprintf(fileID, 'ndof,%d\n', ndof);
fprintf(fileID, 'nnode,%d\n', nnode);

%% node printing

x = linspace(xL, xR, n_elem + 1);
y = linspace(yL, yR, n_elem + 1);

[X, Y] = meshgrid(x, y);
Z = zeros(n_elem + 1, n_elem + 1);

surf(X, Y, Z);

node_coords = zeros(nnode, 2);

for i = 0:n_elem
    x_coords = X(i + 1, :)';
    y_coords = Y(i + 1, :)';
    ki = 1;

    for j = (((n_elem + 1) * i) + 1):((n_elem + 1) * (i + 1))
        node_coords(j, 1) = x_coords(ki);
        node_coords(j, 2) = y_coords(i + 1);
        ki = ki + 1;
        fprintf(fileID, '%d,%f,%f\n', j, node_coords(j, 1), node_coords(j, 2));
    end

end

%% matdata

fprintf(fileID, 'matdata,%d\n', matdata);

for i = 1:matdata
    fprintf(fileID, '%d,%.8f,%f,%f,%f,%f\n', i, k, ux, uy, s, f);
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
conn_table = zeros(n_elem + 1, n_elem + 1);
xx = n_elem + 1;

for i = 0:n_elem

    for j = 1:(n_elem + 1)
        conn_table(i + 1, j) = xx * i + j;
    end

end

% print element connectivity table
xx = 1;

for i = 0:n_elem - 1
    l1 = conn_table(i + 1, :)';
    l2 = conn_table(i + 2, :)';

    for j = 0:n_elem - 1
        fprintf(fileID, '%d,%d,%d,%d,%d,%d,%d\n', xx, matdata, elemdata, l1(j + 1), l1(j + 2), l2(j + 2), l2(j + 1));
        xx = xx + 1;
    end

end

%% boundary layers

fprintf(fileID, "nbound,%d\n", nbound);

% kk = n_elem + 1;
% rev_conn_table = zeros(n_elem + 1, n_elem + 1);
% for i = 0:n_elem
% 	for j = 1:(n_elem + 1)
% 		rev_conn_table(kk, j) = conn_table(i+1, j);
% 	end
% 	kk = kk - 1;
% end

% for x = 0, left edge
bc_ycoords = Y(:, 1); % all values in first column

idx_y = find(bc_ycoords == 0.7); % find the index for y values starting from 0.7
% print values for x (0, 0.7) to (0, 1)
node_nums = conn_table(:, 1);
b_val = 1.0; % boundary value

for i = 1:length(bc_ycoords)

    if (i >= idx_y)
        fprintf(fileID, "%d,%d,%f\n", node_nums(i), elemdata, b_val);
    else
        fprintf(fileID, "%d,%d,%f\n", node_nums(i), elemdata, 0.0);
    end

end

% for y = 0, top edge
bc_xcoords = X(end, 2:end); % final row, 2nd col to end col because 0,1 is already done

% print values for x (0, 1) to (1, 1)
node_nums = conn_table(end, 2:end);
b_val = 1.0; % boundary value

for i = 1:length(bc_xcoords)
    fprintf(fileID, "%d,%d,%f\n", node_nums(i), elemdata, b_val);
end

% bottom edge
bc_xcoords = X(1, 2:end); % final row, 2nd col to end col because 0,1 is already done

% print values for x (0, 0) to (1, 0)
node_nums = conn_table(1, 2:end);
b_val = 0.0; % boundary value

for i = 1:length(bc_xcoords)
    fprintf(fileID, "%d,%d,%f\n", node_nums(i), elemdata, b_val);
end

% top edge
bc_ycoords = Y((2:(end - 1)), end); % final row, 2nd col to end col because 0,1 is already done

% print values for x (1, 0) to (1, 1)
node_nums = conn_table((2:(end - 1)), end);
b_val = 0.0; % boundary value

for i = 1:length(bc_ycoords)
    fprintf(fileID, "%d,%d,%f\n", node_nums(i), elemdata, b_val);
end

%% forcing term
fprintf(fileID, "nforce,%d\n", nforce);

% end
