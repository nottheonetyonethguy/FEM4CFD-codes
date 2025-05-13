%% 1D steady state advection diffusion
clear;

xL = 0;
xR = 1;
mu = 0.0080;
c = 2;
f = 0;
nelem = 10;

L = xR - xL;
he = L / nelem;

% boundary conditions
uL = 0;
uR = 1;

Pe = (c * he) / (2 * mu);

nGP = 2;
[gpts, gwts] = get_Gausspoints_1D(nGP);
alpha = 1 / tanh(Pe) - 1 / Pe;
tau = (he / (2 * c)) * alpha;

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
soln_full_ppv = zeros(totaldof, 1);
soln_full1 = zeros(totaldof, 1);
%% Solution

Ke_diffusion = (mu / he) * [1, -1; -1, 1];
Ke_advection = (c / 2) * [-1, 1; -1, 1];

for iter = 1:9

    Kglobal_g = zeros(totaldof, totaldof);
    Kglobal_pg = zeros(totaldof, totaldof);
    Kglobal_supg = zeros(totaldof, totaldof);
    Kglobal_ppv = zeros(totaldof, totaldof);
    Fglobal_g = zeros(totaldof, 1);
    Fglobal_pg = zeros(totaldof, 1);
    Fglobal_supg = zeros(totaldof, 1);
    Fglobal_ppv = zeros(totaldof, 1);

    Kglobal1 = zeros(totaldof, totaldof);
    Fglobal1 = zeros(totaldof, 1);

    for elnum = 1:nelem
        elem_dofs = elem_dof_conn(elnum, :);
        Klocal = zeros(2, 2);
        Klocal = galerkinApproximation(c, mu, he, nGP, gpts, gwts, Klocal);
	    Kglobal_g(elem_dofs, elem_dofs) = Kglobal_g(elem_dofs, elem_dofs) + Klocal; % galerkin approximation
	    Klocal_pg = Klocal;
	    Klocal_pg = petrovGalerkin(c, mu, he, alpha, nGP, gpts, gwts, Klocal_pg);
	    Kglobal_pg(elem_dofs, elem_dofs) = Kglobal_pg(elem_dofs, elem_dofs) + Klocal_pg;
	    Klocal = supg(c, mu, he, alpha, tau, nGP, gpts, gwts, Klocal);
	    Kglobal_supg(elem_dofs, elem_dofs) = Kglobal_supg(elem_dofs, elem_dofs) + Klocal;
	    
	    Klocal = ppv(c, mu, he, alpha, tau, 0, nGP, gpts, gwts, elem_dofs, node_coords, Klocal);
	    Kglobal_ppv(elem_dofs, elem_dofs) = Kglobal_ppv(elem_dofs, elem_dofs) + Klocal;
        % profs solution
	    [Klocal1, Flocal1] = calcStiffnessAndForce_1D2noded_AdvectionDiffusionReaction(elem_dofs, node_coords, c, mu, 1, f, soln_full, 1, 1.0);
	    Kglobal1(elem_dofs, elem_dofs) = Kglobal1(elem_dofs, elem_dofs) + Klocal1;
	    Fglobal1(elem_dofs, 1) = Fglobal1(elem_dofs, 1) + Flocal1;
    end

    Fglobal_g = forceVector(Kglobal_g, Fglobal_g, iter, uL, uR, totaldof);
    Fglobal_pg = forceVector(Kglobal_pg, Fglobal_pg, iter, uL, uR, totaldof);
    Fglobal_supg = forceVector(Kglobal_supg, Fglobal_supg, iter, uL, uR, totaldof);
    Fglobal_ppv = forceVector(Kglobal_ppv, Fglobal_ppv, iter, uL, uR, totaldof);
    Fglobal1 = forceVector(Kglobal1, Fglobal1, iter, uL, uR, totaldof);

    rNorm = norm(Fglobal_g);

    if (rNorm < 1.0e-10)
	break;
    end

    Kglobal_g = stiffnessMatrix(Kglobal_g, totaldof);
    Kglobal_pg = stiffnessMatrix(Kglobal_pg, totaldof);
    Kglobal_supg = stiffnessMatrix(Kglobal_supg, totaldof);
    Kglobal_ppv = stiffnessMatrix(Kglobal_ppv, totaldof);
    Kglobal1 = stiffnessMatrix(Kglobal1, totaldof);

    soln_incr = Kglobal_g \ Fglobal_g;
    soln_incr_pg = Kglobal_pg \ Fglobal_pg;
    soln_incr_supg = Kglobal_supg \ Fglobal_supg;
    soln_incr_ppv = Kglobal_ppv \ Fglobal_ppv;
    soln_incr1 = Kglobal1 \ Fglobal1;

    soln_full = soln_full + soln_incr;
    soln_full_pg = soln_full_pg + soln_incr_pg;
    soln_full_supg = soln_full_supg + soln_incr_supg;
    soln_full_ppv = soln_full_ppv + soln_incr_ppv;
    soln_full1 = soln_full1 + soln_incr1;

end

u_analytical = analyticalSolution(node_coords, c, mu, L);

f1 = figure;
f2 = figure;
figure(f1);
plot(node_coords, soln_full, 'ro-', 'DisplayName', 'Galerkin');
hold on;
plot(node_coords, soln_full_pg, 'k--', 'DisplayName', 'Petrov-Galerkin');
plot(node_coords, soln_full_supg, 'b *-', 'DisplayName', 'SUPG');
plot(node_coords, soln_full_ppv, 'g--', 'DisplayName', 'Discontinuity Correction');
plot(node_coords, u_analytical', 'k-', 'LineWidth', 2, 'DisplayName', 'Analytical Soliution')
xlabel("Node Coordinates")
ylabel("Values")
title("1D steady state advection diffusion equation")
legend('Location', 'best');

figure(f2);
plot(node_coords, soln_full1, 'ko-', 'DisplayName', 'Function Solution');
hold on;
plot(node_coords, soln_full_ppv, 'r--', 'DisplayName', 'Discontinuity Correction');
xlabel("Node Coordinates")
ylabel("Values")
title("1D Steady State Advection with functions")
legend('Location', 'best');

function Klocal_g = galerkinApproximation(a, mu, h, nGP, gpts, gwts, Klocal)
    Klocal_g = zeros(2, 2);
    for gp = 1:nGP
	xi = gpts(gp);
	wt = gwts(gp);
	N = [0.5 * (1 - xi), 0.5 * (1 + xi)];
	dNdxi = [-0.5, 0.5];
	Jac = h / 2;
	dNdx = dNdxi / Jac;
	% advection
	Klocal = Klocal + (a * N' * dNdx) * Jac * wt;
	%diffusion
	Klocal = Klocal + (mu * dNdx' * dNdx) * Jac * wt;
    end
    Klocal_g = Klocal_g + Klocal;
end

function Klocal_pg = petrovGalerkin(a, mu, h, alpha, nGP, gpts, gwts, Klocal)
    Klocal_pg = zeros(2, 2);
    for gp = 1:nGP
	    xi = gpts(gp);
	    wt = gwts(gp);

	    N = [0.5 * (1 - xi), 0.5 * (1 + xi)];
	    dNdxi = [-0.5, 0.5];

	    Jac = h / 2;

	    dNdx = dNdxi / Jac;

	    % petrov galerkin advection
	    Klocal = Klocal + (alpha * a * dNdx' * dNdx) * Jac * wt;
    end
    Klocal_pg = Klocal_pg + Klocal;
end

function Klocal_supg = supg(a, mu, h, alpha, tau, nGP, gpts, gwts, Klocal)
    Klocal_supg = zeros(2, 2);

	for gp = 1:nGP
	    xi = gpts(gp);
	    wt = gwts(gp);

	    N = [0.5 * (1 - xi), 0.5 * (1 + xi)];
	    dNdxi = [-0.5, 0.5];

	    Jac = h / 2;

	    dNdx = dNdxi / Jac;

	    % stabilization
	    Klocal = Klocal + tau * a ^ 2 * (dNdx' * dNdx) * Jac * wt;
	end

	Klocal_supg = Klocal_supg + Klocal;

end

function Klocal_ppv = ppv(a, mu, h, alpha, tau, f, nGP, gpts, gwts, elem_dofs, soln_full_ppv, Klocal)
    Klocal_ppv = zeros(2, 2);
    n1 = elem_dofs(1);
    n2 = elem_dofs(2);
    
    u1 = soln_full_ppv(n1);
    u2 = soln_full_ppv(n2);
    u = [u1 u2];

	for gp = 1:nGP
	    xi = gpts(gp);
	    wt = gwts(gp);

	    N = [0.5 * (1 - xi), 0.5 * (1 + xi)];
	    dNdxi = [-0.5, 0.5];

	    Jac = h / 2;

	    dNdx = dNdxi / Jac;
	    
	    du = dNdx * u'; 

	    if f == 0
		res_ratio = abs(a);
	    else
		res_ratio = (abs(a * du - f) / abs(du));
	    end

	    chi = 2 / (abs(f) * h + 2 * abs(a));

	    kadd = max(0, ((abs(a - tau * a * f + tau * a * abs(f)) * h / 2) - (mu + tau * a * a) + (f + tau * f * abs(f) * h * h / 6)));
	    
	    if a > 0.0 
		res_ratio = -res_ratio;
	    end

	    % stabilization
	    Klocal = Klocal + chi * res_ratio * kadd * (dNdx' * dNdx) * Jac * wt;
	end

	Klocal_ppv = Klocal_ppv + Klocal;

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

function u_analytical = analyticalSolution(x, a, mu, L)
    Pe = (a * L) / mu;

    u_analytical = (exp(a * x * mu) - 1) / (exp(Pe) - 1);

end
