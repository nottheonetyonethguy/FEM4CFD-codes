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
u = init_soln;
plot_results(x, u, gamma, elem_node_conn, node_coords);
u_l = u(:, 1); u_r = u(:, end);

soln_initial = soln_applied;
soln_full_prev = soln_initial;

soln_applied = soln_applied(:);

for loadstep = 1:maxsteps
	
    u_p = soln_applied;
	Kglobal = zeros(4 * nnode, 4 * nnode);
	Mglobal = zeros(4 * nnode, 4 * nnode);
	
	for elnum = 1:nelem
		node_idx = elem_node_conn(elnum, :);
		U = u(:, node_idx);
		% caclulate the local stiffness matrix
		[Klocal, Mlocal, K_supg, M_supg] = stiffnessAndForce_quad4_euler2D(elem_node_conn(elnum, :), node_coords, elemData, matData, U, dt, soln_initial, soln_full_prev);
		
		% get the indices for assembly
		inds4Assy = elem_dof_conn(elnum, :);
		
		% assemble
		Kglobal(inds4Assy, inds4Assy) = Kglobal(inds4Assy, inds4Assy) + Klocal + K_supg;
		Mglobal(inds4Assy, inds4Assy) = Mglobal(inds4Assy, inds4Assy) + Mlocal + M_supg;
    end
    
    LHS = (1/dt) * Mglobal + Kglobal;
	RHS = (1/dt) * Mglobal * u(:);
    
	for i = 1:max(size(dofs_fixed))
		dof = dofs_fixed(i);
		RHS = RHS - LHS(:, dof) * soln_applied(dof);
	end
	u_new = LHS(dofs_free, dofs_free) \ RHS(dofs_free);

    for i = 1:max(size(dofs_free))
        dof = dofs_free(i);
        u_p(dof) = u_p(dof) + u_new(i);
    end
    u = reshape(u_p, ndof, nnode);

	% fileName = sprintf("%s%d%s", "advdiff2D-", loadstep, "-.vtk");
	% writeoutputvtk(ndim, nelem, nnode, npelem, ndof, node_coords, elem_node_conn, soln_initial, fileName);
end

plot_results(x, u, gamma, elem_node_conn, node_coords);

