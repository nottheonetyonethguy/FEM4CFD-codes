clear;
clc;

% Step 1: Pre-prcessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fname = "QuadMesh_mesh2_Structured.inp";
% fname = "data0.inp";
% fname = "data30.inp";
% fname = "data60.inp";
% fname = "data45.inp";
fname = "data90.inp";
% fname = "data_p2.inp";

[ndim, ndof, nnode, nelem, npelem, totaldof, node_coords, elem_node_conn, matData, elemData, elem_dof_conn, dofs_fixed, dofs_free, soln_applied, soln_full, force_applied] = processfile_FEM(fname);

% Step 2 - Solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = 0.1;

maxsteps = 1.0 / dt;

load_factor = 0.0;

soln_full = soln_full * 0.0;
soln_full_prev = soln_full;

for loadstep = 1:maxsteps
    disp(sprintf("\nLoad step %d \n", loadstep));
    disp(sprintf("================================\n"));
    load_factor = load_factor + dt;

    soln_applied_step = dt * soln_applied;

    for iter = 1:4

        Kglobal = zeros(totaldof, totaldof);
        Fglobal = zeros(totaldof, 1);

        for elnum = 1:nelem
            % caclulate the local stiffness matrix
            [Klocal, Flocal] = stiffnessAndForce_quad4_AdvDiff2D(elem_node_conn(elnum, :), node_coords, force_applied, elemData, matData, soln_full, soln_full_prev);

            % get the indices for assembly
            inds4Assy = elem_dof_conn(elnum, :);

            % Assemble
            Kglobal(inds4Assy, inds4Assy) = Kglobal(inds4Assy, inds4Assy) + Klocal;
            Fglobal(inds4Assy, 1) = Fglobal(inds4Assy, 1) + Flocal;
        end

        % apply boundary conditions
        if (iter == 1)

            for i = 1:max(size(dofs_fixed))
                dof = dofs_fixed(i);
                Fglobal = Fglobal - Kglobal(:, dof) * soln_applied_step(dof);
            end

            for i = 1:max(size(dofs_fixed))
                dof = dofs_fixed(i);
                Fglobal(dof, 1) = soln_applied_step(dof);
            end

        else

            for i = 1:max(size(dofs_fixed))
                Fglobal(dofs_fixed(i), 1) = 0.0;
            end

        end

        rNorm = norm(Fglobal);

        disp(sprintf("Iteration %d ... norm = %E \n", iter, rNorm));

        if (rNorm < 1.0e-10)
            break;
        end

        % Replace the rows and colums corresponding to the specified DOFs with zeros
        % and add a diagonal entry of 1
        for i = 1:max(size(dofs_fixed))
            dof = dofs_fixed(i);

            Kglobal(dof, :) = zeros(totaldof, 1);
            Kglobal(:, dof) = zeros(totaldof, 1);
            Kglobal(dof, dof) = 1.0;
        end

        % solve for unknown DOFs
        soln_incr = Kglobal \ Fglobal;

        % populate the full solution array
        soln_full = soln_full + soln_incr;
    end

    soln_full_prev = soln_full;

    fileName = sprintf("%s%d%s", "advdiff2D-", loadstep, "-.vtk");
    writeoutputvtk(ndim, nelem, nnode, npelem, ndof, node_coords, elem_node_conn, soln_full, fileName);
end

% Step 3 - Post-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%trimesh(elem_node_conn, node_coords(:,1), node_coords(:,2))

%figure(2)
%trimesh(elem_node_conn, node_coords(:,1), node_coords(:,2), soln_full)
trisurf(elem_node_conn, node_coords(:, 1), node_coords(:, 2), soln_full, 'edgecolor', 'k', 'facecolor', 'interp');
colormap jet;
view(2); colorbar;
