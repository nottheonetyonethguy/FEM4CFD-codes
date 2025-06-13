%% 1D steady state advection diffusion
clear;

xL = 0;
xR = 1;
mu = 0.0016;
c = 2;
f = 0;
nelem = 50;

L = xR - xL;
he = L / nelem;

% boundary conditions
uL = 0;
uR = 1;

Pe = (c * he) / (2 * mu);

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
soln_full_pg = zeros(totaldof, 1);
soln_full_supg = zeros(totaldof, 1);
soln_full1 = zeros(totaldof, 1);
%% Solution

Ke_diffusion = (mu / he) * [1, -1; -1, 1];
Ke_advection = (c / 2) * [-1, 1; -1, 1];

for iter = 1:9
		
	Kglobal_g = zeros(totaldof, totaldof);
	Kglobal_pg = zeros(totaldof, totaldof);
	Kglobal_supg = zeros(totaldof, totaldof);
	Fglobal_g = zeros(totaldof, 1);
	Fglobal_pg = zeros(totaldof, 1);
	Fglobal_supg = zeros(totaldof, 1);
	
	Kglobal1 = zeros(totaldof, totaldof);
	Fglobal1 = zeros(totaldof, 1);
	
	Kglobal_g = galerkinApproximation(c, mu, he, elem_dof_conn, totaldof, nelem);
	
	Kglobal_pg = petrovGalerkin(c, mu, he, Pe, elem_dof_conn, totaldof, nelem);
	
	% Kglobal_supg = supg(c, mu, he, Pe, elem_dof_conn, totaldof, nelem);
	
	% for elnum = 1:nelem
	%     % element connectivity
	%     elem_dofs = elem_dof_conn(elnum, :);
	%     % profs solution
	%     [Klocal1, Flocal1] = calcStiffnessAndForce_1D2noded_AdvectionDiffusionReaction(elem_dofs, node_coords, c, mu, 1, f, soln_full, 1, 1.0);
	%     Kglobal1(elem_dofs, elem_dofs) = Kglobal1(elem_dofs, elem_dofs) + Klocal1;
	%     Fglobal1(elem_dofs, 1) = Fglobal1(elem_dofs, 1) + Flocal1;
	% end
	
	Fglobal_g = forceVector(Kglobal_g, Fglobal_g, iter, uL, uR, totaldof);
	Fglobal_pg = forceVector(Kglobal_pg, Fglobal_pg, iter, uL, uR, totaldof);
	% Fglobal_supg = forceVector(Kglobal_supg, Fglobal_supg, iter, uL, uR, totaldof);
	% Fglobal1 = forceVector(Kglobal1, Fglobal1, iter, uL, uR, totaldof);
	
	rNorm = norm(Fglobal_g);
	
	if (rNorm < 1.0e-10)
		break;
	end
	
	Kglobal_g = stiffnessMatrix(Kglobal_g, totaldof);
	Kglobal_pg = stiffnessMatrix(Kglobal_pg, totaldof);
	% Kglobal_supg = stiffnessMatrix(Kglobal_supg, totaldof);
	% Kglobal1 = stiffnessMatrix(Kglobal1, totaldof);
	
	soln_incr = Kglobal_g \ Fglobal_g;
	soln_incr_pg = Kglobal_pg \ Fglobal_pg;
	% soln_incr_supg = Kglobal_supg \ Fglobal_supg;
	% soln_incr1 = Kglobal1 \ Fglobal1;
	
	soln_full = soln_full + soln_incr;
	soln_full_pg = soln_full_pg + soln_incr_pg;
	% soln_full_supg = soln_full_supg + soln_incr_supg;
	% soln_full1 = soln_full1 + soln_incr1;
	
end

u_analytical = analyticalSolution(node_coords, c, mu, he, L);

% subplot(1, 2, 1)
plot(node_coords, soln_full, 'ro-', 'DisplayName', 'Galerkin');
hold on;
plot(node_coords, soln_full_pg, 'k--', 'DisplayName', 'Petrov-Galerkin');
% plot(node_coords, soln_full_supg, 'b *-', 'DisplayName', 'SUPG');
plot(node_coords, u_analytical', 'k-', 'LineWidth', 2, 'DisplayName', 'Analytical Soliution')
xlabel("Node Coordinates")
ylabel("Values")
title("1D steady state advection diffusion equation")
legend('Location', 'best');

% subplot(1, 2, 2)
% plot(node_coords, soln_full1, 'ko-');
% xlabel("Node Coordinates")
% ylabel("Values")
% title("1D Steady State Advection with functions")

function Kglobal_g = galerkinApproximation(a, mu, h, elem_dof_conn, totaldof, nelem)
Kglobal_g = zeros(totaldof, totaldof);
K_advection = (a / 2) * [-1 1; -1 1];
K_diffusion = (mu / h) * [1 -1; -1 1];

for elnum = 1:nelem
	elem_dofs = elem_dof_conn(elnum, :);
	Klocal = Kglobal_g(elem_dofs, elem_dofs);
	Klocal = Klocal + K_advection + K_diffusion;
	Kglobal_g(elem_dofs, elem_dofs) = Klocal;
end

end

function Kglobal_pg = petrovGalerkin(a, mu, h, Pe, elem_dof_conn, totaldof, nelem)
Kglobal_pg = zeros(totaldof, totaldof);
nGP = 2;
[gpts, gwts] = get_Gausspoints_1D(nGP);
% alpha = 1;
alpha = 1 / tanh(Pe) - 1 / Pe;

for elnum = 1:nelem
	elem_dofs = elem_dof_conn(elnum, :);
	Klocal = zeros(2, 2);
	
	for gp = 1:nGP
		xi = gpts(gp);
		wt = gwts(gp);
		
		N = [0.5 * (1 - xi), 0.5 * (1 + xi)];
		dNdxi = [-0.5, 0.5];
		
		Jac = h / 2;
		
		dNdx = dNdxi / Jac;
		
		% advection
		Klocal = Klocal + (a * N' * dNdx) * Jac * wt;
		Klocal = Klocal + (alpha * a * dNdx' * dNdx) * Jac * wt;
		
		% diffusion
		Klocal = Klocal + mu * (dNdx' * dNdx) * Jac * wt;
	end
	
	Kglobal_pg(elem_dofs, elem_dofs) = Kglobal_pg(elem_dofs, elem_dofs) + Klocal;
end

end

function Kglobal_supg = supg(a, mu, h, Pe, elem_dof_conn, totaldof, nelem)
Kglobal_supg = zeros(totaldof, totaldof);

nGP = 2;
[gpts, gwts] = get_Gausspoints_1D(nGP);

alpha = 1 / tanh(Pe) - 1 / Pe;
tau = (h / (2 * a)) * alpha;

for elnum = 1:nelem
	elem_dofs = elem_dof_conn(elnum, :);
	Klocal = zeros(2, 2);
	
	for gp = 1:nGP
		xi = gpts(gp);
		wt = gwts(gp);
		
		N = [0.5 * (1 - xi), 0.5 * (1 + xi)];
		dNdxi = [-0.5, 0.5];
		
		Jac = h / 2;
		
		dNdx = dNdxi / Jac;
		
		% advection
		Klocal = Klocal + (a * N' * dNdx) * Jac * wt;
		
		% diffusion
		Klocal = Klocal + mu * (dNdx' * dNdx) * Jac * wt;
		
		% stabilization
		Klocal = Klocal + tau * a ^ 2 * (dNdx' * dNdx) * Jac * wt;
	end
	
	Kglobal_supg(elem_dofs, elem_dofs) = Kglobal_supg(elem_dofs, elem_dofs) + Klocal;
end

end

function Kglobal = stiffnessMatrix(Kglobal, totaldof)
Kglobal(1, :) = zeros(totaldof, 1);
Kglobal(:, 1) = zeros(totaldof, 1);
Kglobal(1, 1) = 1.0;

Kglobal(end, :) = zeros(totaldof, 1);
Kglobal(:, end) = zeros(totaldof, 1);
Kglobal(end, end) = 1.0;
end

function Fglobal = forceVector(Kglobal, Fglobal, iter, uL, uR, totaldof)

if iter == 1
	Fglobal = Fglobal - Kglobal(:, 1) * uL;
	Fglobal = Fglobal - Kglobal(:, totaldof) * uR;
	Fglobal(1, 1) = uL;
	Fglobal(end, 1) = uR;
else
	Fglobal(1, 1) = 0.0;
	Fglobal(end, 1) = 0.0;
	
end

end

function u_analytical = analyticalSolution(x, a, mu, he, L)
Pe = (a * he) / (2 * mu);
k = mu;

u_analytical = (exp(a * x * mu) - 1) / (exp(Pe) - 1);

end
