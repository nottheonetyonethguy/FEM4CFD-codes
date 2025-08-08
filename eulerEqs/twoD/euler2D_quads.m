clear;
clc;
close all;

% Step 1: Pre-prcessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname = "data0.inp";

[ndim, ndof, nnode, nelem, npelem, totaldof, node_coords, elem_node_conn, matData, elemData, elem_dof_conn, init_soln, dofs_fixed, dofs_free, soln_applied, soln_initial] = processfile_FEM(fname);

% Step 2 - Solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_l = 0; x_r = 1;
x = linspace(x_l, x_r, 21);
y = linspace(0, 0.1, 4);
gamma = matData(1, 1);

dt = 0.001;

maxsteps = 0.2 / dt;
maxsteps = 2;
u = init_soln(:, 1:21);
plot_results(x, u, gamma);
trisurf(elem_node_conn, x, y, )
u_l = u(:, 1); u_r = u(:, end);

soln_initial = soln_applied;
soln_full_prev = soln_initial;

for loadstep = 1:maxsteps
		
		Kglobal = zeros(4 * nnode, 4 * nnode);
		Mglobal = zeros(4 * nnode, 4 * nnode);
		
		for elnum = 1:nelem
				node_idx = elem_node_conn(elnum, :);
				U = u(:, node_idx);
				% caclulate the local stiffness matrix
				[Klocal, Mlocal] = stiffnessAndForce_quad4_euler2D(elem_node_conn(elnum, :), node_coords, elemData, matData, U, dt, soln_initial, soln_full_prev);
				
				% get the indices for assembly
				inds4Assy = elem_dof_conn(elnum, :);
				
				% Assemble
				Kglobal(inds4Assy, inds4Assy) = Kglobal(inds4Assy, inds4Assy) + Klocal;
				Mglobal(inds4Assy, inds4Assy) = Mglobal(inds4Assy, inds4Assy) + Mlocal;
		end
		dofs_free = 5:(nnode * 4 - 4);
		LHS = (1/dt) * Mglobal + Kglobal;
		RHS = (1/dt) * Mglobal * u(:);
		
		% RHS = RHS - LHS(:, 1) * U(1, 1) - LHS(:, 2) * U(2, 1) - LHS(:, 3) * U(3, 1) - LHS(:, 4) * U(4, 1);
		% RHS = RHS - LHS(:, end-3) * U(1, end) - LHS(:, end-2) * U(2, end) - LHS(:, end-1) * U(3, end) - LHS(:, end) * U(4, end);
		
		u_new = LHS(dofs_free, dofs_free) \ RHS(dofs_free);
		u(:, 2:end-1) = reshape(u_new, ndof, nnode-2);
		u(:, 1) = u_l; u(:, end) = u_r;
		
		
		% fileName = sprintf("%s%d%s", "advdiff2D-", loadstep, "-.vtk");
		% writeoutputvtk(ndim, nelem, nnode, npelem, ndof, node_coords, elem_node_conn, soln_initial, fileName);
end

% Step 3 - Post-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%trimesh(elem_node_conn, node_coords(:,1), node_coords(:,2))

%figure(2)
%trimesh(elem_node_conn, node_coords(:,1), node_coords(:,2), soln_full)
% trisurf(elem_node_conn, node_coords(:, 1), node_coords(:, 2), soln_initial, 'edgecolor', 'k', 'facecolor', 'interp');
% colormap jet;
% view(2); colorbar;

plot_results(x, u, gamma);

