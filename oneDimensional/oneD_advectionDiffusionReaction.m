%% 1D steady state advection diffusion
clc; clear;

xL = 0;
xR = 1;
mu = 1 % for Pe = 1;
c = 2;

a = 1;
nelem = 9;

L = xR - xL;
he = L / nelem;

% boundary conditions
uL = sin(-a);
uR = sin(a);

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
%% Solution

Ke_diffusion = (mu / he) * [1, -1; -1, 1];
Ke_advection = (c / 2) * [-1, 1; -1, 1];

for iter = 1:9
    Kglobal = zeros(totaldof, totaldof);
    Fglobal = zeros(totaldof, 1);

    for elnum = 1:nelem
        % element connectivity
        elem_dofs = elem_dof_conn(elnum, :);
        % my solution
        Klocal = Kglobal(elem_dofs, elem_dofs);
        Klocal = Klocal + Ke_diffusion + Ke_advection;
        Kglobal(elem_dofs, elem_dofs) = Klocal;

        % for sine force term
        xa = node_coords(elnum);
        xb = node_coords(elnum + 1);

        Flocal = zeros(2, 1);
        Flocal(1) = ((xb - xa) / he) * (c * sin(a * xa) + mu * a * cos(a * xa));
        Flocal(2) = ((xb - xa) / he) * (c * sin(a * xb) - mu * a * cos(a * xb));

        Kglobal(elem_dofs, elem_dofs) = Kglobal(elem_dofs, elem_dofs) + Klocal;
        Fglobal(elem_dofs, 1) = Fglobal(elem_dofs, 1) + Flocal;
    end

    if iter == 1
        Fglobal = Fglobal - Kglobal(:, 1) * uL;
        Fglobal = Fglobal - Kglobal(:, totaldof) * uR;
        Fglobal(1, 1) = uL;
        Fglobal(end, 1) = uR;
    else
        Fglobal(1, 1) = 0.0;
        Fglobal(end, 1) = 0.0;
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

    soln_incr = Kglobal \ Fglobal;
    soln_full = soln_full + soln_incr;
end

plot(node_coords, soln_full, 'ro-');
xlabel("Node Coordinates")
ylabel("Values")
title("1D steady state advection diffusion equation with force term")
