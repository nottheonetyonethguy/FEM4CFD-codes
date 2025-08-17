%% 2D steady state advection diffusion

close all;
clear;

nelem = 10;
x = linspace(0, 1, nelem+1);
y = linspace(0, 1, nelem+1);

[X, Y] = meshgrid(x, y);
Z = zeros(11, 11);

% surf(X, Y, Z);

node_coords = zeros((nelem + 1) ^ 2, 2);

for i=0:nelem
	x_coords = X(i + 1, :)';
	y_coords = Y(i + 1, :)';
	k = 1;
	for j = (((nelem + 1) * i) + 1):((nelem + 1) * (i + 1))
		node_coords(j, 1) = x_coords(k);
		node_coords(j, 2) = y_coords(i + 1);
		k = k + 1;
	end
end

fileID = fopen('test.inp', 'a+');
n_elem = 25;
conn_table = zeros(n_elem + 1, n_elem + 1);
xx = 26;
for i = 0:n_elem
		for j = 1:(n_elem + 1)
				conn_table(i+1, j) = xx * i + j;
		end
end

matdata = 1;
elemdata = 1;
% print element connectivity table
xx = 1;
for  i=0:n_elem-1
		l1 = conn_table(i+1, :)';
		l2 = conn_table(i+2, :)';
		for j = 0:n_elem-1
				fprintf(fileID, '%d,%d,%d,%d,%d,%d,%d\n', xx, matdata, elemdata, l1(j+1), l1(j+2), l2(j+2), l2(j+1));
				xx = xx + 1;
		end
end


