%% 1D steady state advection diffusion
clc; clear; close all;


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

elem_node_conn = [1:nelem; 2: nnode]';
elem_dof_conn = elem_node_conn;

dofs_full = 1:totaldof;
dofs_fixed = [1, totaldof];
dofs_free = setdiff(dofs_full, dofs_fixed);

% solution array
soln_full = zeros(totaldof, 1);
soln_full(1, 1) = uL;
soln_full(totaldof, 1) = uR;

%% Solution

Kglobal = zeros(totaldof, totaldof);
Fglobal = zeros(totaldof, 1);

Ke_diffusion = (mu / he) * [1, -1; -1, 1];
Ke_advection = (c / 2) * [-1, 1; -1, 1];

for elnum = 1:nelem
    % element connectivity
    elem_dofs = elem_dof_conn(elnum, :);
    Klocal = Kglobal(elem_dofs, elem_dofs);
    Klocal = Klocal + Ke_diffusion + Ke_advection;
    Kglobal(elem_dofs, elem_dofs) = Klocal;
    Fglobal(elem_dofs,1) = 0; 
end

% left bc
Fglobal = Fglobal - Kglobal(:, 1) * uL;

% right bc
Fglobal = Fglobal - Kglobal(:, totaldof) * uR;

% remove fixed dofs
Kglobal = Kglobal(dofs_free, dofs_free);
Fglobal = Fglobal(dofs_free, 1);

soln_free_dofs = Kglobal \ Fglobal;

soln_full(dofs_free) = soln_free_dofs;

plot(node_coords, soln_full, 'ro-');

