%% 1D steady state advection diffusion
clc; clear;

xL = 0;
xR = 1;
mu = 1;
c = 2;
f = 0;
nelem = 9;

L = xR - xL;
he = L / nelem;

% boundary conditions
uL = 1;
uR = 0;

nnode = nelem + 1;

ndof = 1;

totaldof = nnode * ndof;

node_coords = linspace(xL, xR, nnode);

elem_node_conn = [1:nelem; 2:nnode]';
elem_dof_conn = elem_node_conn;

dofs_full = 1:totaldof;
dofs_fixed = [1, totaldof];
dofs_free = setdiff(dofs_full, dofs_fixed);

% solution array
soln_full = zeros(totaldof, 1);
soln_full1 = zeros(totaldof, 1);
%% Solution

Kglobal = zeros(totaldof, totaldof);
Fglobal = zeros(totaldof, 1);

Kglobal1 = zeros(totaldof, totaldof);
Fglobal1 = zeros(totaldof, 1);

Ke_diffusion = (mu / he) * [1, -1; -1, 1];
Ke_advection = (c / 2) * [-1, 1; -1, 1];

for iter = 1:9

    for elnum = 1:nelem
        % element connectivity
        elem_dofs = elem_dof_conn(elnum, :);
        % profs solution
        [Klocal1, Flocal1] = calcStiffnessAndForce_1D2noded_AdvectionDiffusionReaction(elem_dofs, node_coords, c, mu, 1, f, soln_full, 1, 1.0);
        Kglobal1(elem_dofs, elem_dofs) = Kglobal1(elem_dofs, elem_dofs) + Klocal1;
        Fglobal1(elem_dofs, 1) = Fglobal1(elem_dofs, 1) + Flocal1;
        % my solution
        Klocal = Kglobal(elem_dofs, elem_dofs);
        Klocal = Klocal + Ke_diffusion + Ke_advection;
        Kglobal(elem_dofs, elem_dofs) = Klocal;
    end

    if iter == 1
        Fglobal = Fglobal - Kglobal(:, 1) * uL;
        Fglobal = Fglobal - Kglobal(:, totaldof) * uR;
        Fglobal(1, 1) = uL;
        Fglobal(end, 1) = uR;

        Fglobal1 = Fglobal1 - Kglobal1(:, 1) * uL;
        Fglobal1 = Fglobal1 - Kglobal1(:, totaldof) * uR;
        Fglobal1(1, 1) = uL;
        Fglobal1(end, 1) = uR;
    else
        Fglobal(1, 1) = 0.0;
        Fglobal(end, 1) = 0.0;

        Fglobal1(1, 1) = 0.0;
        Fglobal1(end, 1) = 0.0;
    end

    rNorm = norm(Fglobal);

    if (rNorm < 1.0e-10)
        break;
    end

    Kglobal(1, :) = zeros(totaldof, 1);
    Kglobal(:, 1) = zeros(totaldof, 1);
    Kglobal(1, 1) = 1.0;

    Kglobal(end, :) = zeros(totaldof, 1);
    Kglobal(:, end) = zeros(totaldof, 1);
    Kglobal(end, end) = 1.0;

    Kglobal1(1, :) = zeros(totaldof, 1);
    Kglobal1(:, 1) = zeros(totaldof, 1);
    Kglobal1(1, 1) = 1.0;

    Kglobal1(end, :) = zeros(totaldof, 1);
    Kglobal1(:, end) = zeros(totaldof, 1);
    Kglobal1(end, end) = 1.0;

    soln_incr = Kglobal \ Fglobal;
    soln_incr1 = Kglobal1 \ Fglobal1;

    soln_full = soln_full + soln_incr;
    soln_full1 = soln_full1 + soln_incr;

end

subplot(1, 2, 1)
plot(node_coords, soln_full, 'ro-');

% subplot(1,2,2)
% plot(node_coords, soln_full1, 'ko-');
