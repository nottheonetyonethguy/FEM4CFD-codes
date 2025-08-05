
function [ndim, ndof, nnode, nelem, npelem, totaldof, node_coords, elem_node_conn, matData, elemData, elem_dof_conn, init_soln, dofs_fixed, dofs_free, soln_applied, soln_full, force_applied] = processfile_FEM(fname)
clc
% node_coords: global coordinates of the nodes, x, y, and z
% nelem: total number of elements
% nnode: total number of nodes
% ndof: number of degrees of per node
% elem_node_conn: element-to-node connectivity
% elem_dof_conn: element-to-dof connectivity

fid=fopen(fname,'r');

% ndim

line=fgets(fid);
linestr = strsplit(line, ",");
ndim = int32(str2num(linestr{1,2}))

% ndof

line=fgets(fid);
linestr = strsplit(line, ",");
ndof = int32(str2num(linestr{1,2}))

% nodes

line=fgets(fid);
linestr = strsplit(line, ",");
nnode = int32(str2num(linestr{1,2}))

totaldof = nnode*ndof;

node_coords = zeros(nnode,ndim);
for i=1:nnode
		line = fgets(fid);
		linestr = strsplit(line, ",");
		
		node_coords(i,1) = double(str2num(linestr{1,2}));
		node_coords(i,2) = double(str2num(linestr{1,3}));
		if(ndim == 3)
				node_coords(i,3) = double(str2num(linestr{1,4}));
		end
end


% material data

line=fgets(fid);
linestr = strsplit(line, ",");
nmatData = int32(str2num(linestr{1,2}));

matData = zeros(nmatData, 20);
for i=1:nmatData
		line = fgets(fid);
		linestr = strsplit(line, ",");
		
		for j=1:size(linestr,2)-1
				matData(i,j) = double(str2num(linestr{1,j+1}));
		end
end


% element data

line=fgets(fid);
linestr = strsplit(line, ",");
nelemData = int32(str2num(linestr{1,2}));

elemData = zeros(nelemData, 20);
for i=1:nelemData
		line = fgets(fid);
		linestr = strsplit(line, ",");
		
		for j=1:size(linestr,2)-1
				elemData(i,j) = double(str2num(linestr{1,j+1}));
		end
end


% npelem

line=fgets(fid);
linestr = strsplit(line, ",");
npelem = int32(str2num(linestr{1,2}));


% elements

line=fgets(fid);
linestr = strsplit(line, ",");
nelem = int32(str2num(linestr{1,2}));

elem_node_conn = zeros(nelem, npelem, "int32");
for i=1:nelem
		line = fgets(fid);
		linestr = strsplit(line, ",");
		
		for j=1:npelem
				elem_node_conn(i,j) = int32(str2num(linestr{1,3+j}));
		end
end

% initial conditions of u

line=fgets(fid);
linestr = strsplit(line, ",");
nINIT    = int32(str2num(linestr{1,2}));

gamma = matData(1,1);

ninit = zeros(nINIT, 1, "int32");
init_soln = zeros(4, totaldof);
for i=1:nINIT
		line = fgets(fid);
		linestr = strsplit(line, ",");
		
		n1 = int32(str2num(linestr{1,1}));
		n2 = int32(str2num(linestr{1,2}));
		n3 = double(str2num(linestr{1,3})); % rho
		n4 = double(str2num(linestr{1,4})); rho_u = n3 * n4;% u
		n5 = double(str2num(linestr{1,5})); rho_v = n3 * n5;% v
		n6 = double(str2num(linestr{1,6})); % pressure
		E = (n6 / (gamma - 1)) + 0.5 * n3 * (n4^2 + n5^2);
		init_soln(:, i) = double([n3 rho_u rho_v E]);
end

% Dirichlet boundary conditions

line=fgets(fid);
linestr = strsplit(line, ",");
nDBC    = int32(str2num(linestr{1,2}));

gamma = matData(1,1);

dofs_fixed = zeros(nDBC, 1, "int32");
soln_applied = zeros(4, totaldof);
for i=1:nDBC
	line = fgets(fid);
	linestr = strsplit(line, ",");
	
	n1 = int32(str2num(linestr{1,1}));
	n2 = int32(str2num(linestr{1,2}));
	n3 = double(str2num(linestr{1,3})); % rho
	n4 = double(str2num(linestr{1,4})); rho_u = n3 * n4;% u
	n5 = double(str2num(linestr{1,5})); rho_v = n3 * n5;% v
	n6 = double(str2num(linestr{1,6})); % pressure
	E = (n6 / (gamma - 1)) + 0.5 * n3 * (n4^2 + n5^2);
	dofs_fixed(i) = (n1-1)*ndof+n2;
	soln_applied(:, dofs_fixed(i)) = double([n3 rho_u rho_v E]);
end

dofs_free = setdiff([1:totaldof], dofs_fixed)';
soln_full = soln_applied;


% Force boundary conditions

line=fgets(fid);
linestr = strsplit(line, ",");
nFBC    = int32(str2num(linestr{1,2}));

%dof_force = zeros(nFBC, 1, "int32");
force_applied = zeros(totaldof, 1);
for i=1:nFBC
	line = fgets(fid);
	linestr = strsplit(line, ",");
	
	n1 = int32(str2num(linestr{1,1}));
	n2 = int32(str2num(linestr{1,2}));
	ind = (n1-1)*ndof + n2;
	
	%dof_force(i) = ind;
	force_applied(ind) = double(str2num(linestr{1,3}));
end


% data structures
nsize = npelem*ndof;
elem_dof_conn  = zeros(nelem, nsize);

for e=1:nelem
	count = 1;
	for jj=1:npelem
		ind = ndof*(elem_node_conn(e,jj)-1)+1;
		for kk=1:ndof
			elem_dof_conn(e,count) = ind;
			ind = ind + 1;
			count = count + 1;
		end
	end
	count = count - 1;
end



