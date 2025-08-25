%% 1D steady state advection diffusion
% close all;
clear;

xL = 0;
xR = 1;
c = 1;
f = 0; % changes

nelem = 10;
ZZZZZZZZ
L = xR - xL;
he = L / nelem;

% boundary conditions
uL = 0;
uR = 1;

Pe = 10;
Da = 10; % s * he /c

mu = (c * he) / (2 * Pe);
s = Da * c / he;

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
        Flocal = zeros(2, 1);

        [Klocal, Flocal] = galerkinApproximation(c, mu, he, s, nGP, gpts, gwts, elem_dofs, node_coords, soln_full, Klocal, Flocal);
        Kglobal_g(elem_dofs, elem_dofs) = Kglobal_g(elem_dofs, elem_dofs) + Klocal; % galerkin approximation
        Fglobal_g(elem_dofs, 1) = Fglobal_g(elem_dofs, 1) + Flocal; % galerkin approximation

        [Klocal, Flocal] = supg(c, mu, he, alpha, tau, s, nGP, gpts, gwts, elem_dofs, node_coords, soln_full_supg, Klocal, Flocal);
        Kglobal_supg(elem_dofs, elem_dofs) = Kglobal_supg(elem_dofs, elem_dofs) + Klocal;
        Fglobal_supg(elem_dofs, 1) = Fglobal_supg(elem_dofs, 1) + Flocal; % supg approximation

        [Klocal, Flocal] = ppv(c, mu, he, alpha, tau, s, nGP, gpts, gwts, elem_dofs, node_coords, soln_full_ppv, Klocal, Flocal);
        Kglobal_ppv(elem_dofs, elem_dofs) = Kglobal_ppv(elem_dofs, elem_dofs) + Klocal;
        Fglobal_ppv(elem_dofs, 1) = Fglobal_ppv(elem_dofs, 1) + Flocal; % ppv approximation

        % profs solution
        [Klocal1, Flocal1] = calcStiffnessAndForce_1D2noded_AdvectionDiffusionReaction(elem_dofs, node_coords, c, mu, 1, s, soln_full, 1, 1.0);
        Kglobal1(elem_dofs, elem_dofs) = Kglobal1(elem_dofs, elem_dofs) + Klocal1;
        Fglobal1(elem_dofs, 1) = Fglobal1(elem_dofs, 1) + Flocal1;
    end

    Fglobal_g = forceVector(Kglobal_g, Fglobal_g, iter, uL, uR, totaldof);
    % Fglobal_pg = forceVector(Kglobal_pg, Fglobal_pg, iter, uL, uR, totaldof);
    Fglobal_supg = forceVector(Kglobal_supg, Fglobal_supg, iter, uL, uR, totaldof);
    Fglobal_ppv = forceVector(Kglobal_ppv, Fglobal_ppv, iter, uL, uR, totaldof);
    Fglobal1 = forceVector(Kglobal1, Fglobal1, iter, uL, uR, totaldof);

    rNorm = norm(Fglobal_g);

    if (rNorm < 1.0e-10)
        break;
    end

    Kglobal_g = stiffnessMatrix(Kglobal_g, totaldof);
    Kglobal_supg = stiffnessMatrix(Kglobal_supg, totaldof);
    Kglobal_ppv = stiffnessMatrix(Kglobal_ppv, totaldof);
    Kglobal1 = stiffnessMatrix(Kglobal1, totaldof);

    soln_incr = Kglobal_g \ Fglobal_g;
    soln_incr_supg = Kglobal_supg \ Fglobal_supg;
    soln_incr_ppv = Kglobal_ppv \ Fglobal_ppv;
    soln_incr1 = Kglobal1 \ Fglobal1;

    soln_full = soln_full + soln_incr;
    soln_full_supg = soln_full_supg + soln_incr_supg;
    soln_full_ppv = soln_full_ppv + soln_incr_ppv;
    % soln_full_ppv = max(0, soln_full_ppv);
    soln_full1 = soln_full1 + soln_incr1;

end

u_analytical = analyticalSolution(node_coords, c, mu, s, L, uL, uR);

f1 = figure;
% f2 = figure;
figure(f1);
plot(node_coords, soln_full, 'bo-', 'DisplayName', 'Galerkin');
hold on;
% plot(node_coords, soln_full_pg, 'k--', 'DisplayName', 'Petrov-Galerkin');
plot(node_coords, soln_full_supg, 'r *-', 'DisplayName', 'SUPG');
plot(node_coords, soln_full_ppv, 'k-o', 'DisplayName', 'Discontinuity Correction', 'MarkerFaceColor', 'black');
plot(node_coords, u_analytical', 'k-', 'DisplayName', 'Analytical Soliution')
xlabel("Node Coordinates, x/L")
ylabel("Values, \phi")
title("1D steady state advection diffusion equation")
legend('Location', 'best');

% figure(f2);
% plot(node_coords, soln_full1, 'ko-', 'DisplayName', 'Function Solution');
% hold on;
% plot(node_coords, soln_full_ppv, 'r--', 'DisplayName', 'Discontinuity Correction');
% xlabel("Node Coordinates")
% ylabel("Values")
% title("1D Steady State Advection with functions")
% legend('Location', 'best');

function [Klocal_g, Flocal_g] = galerkinApproximation(a, mu, h, s, nGP, gpts, gwts, elem_dofs, node_coords, soln_full, Klocal, Flocal)
    Klocal_g = zeros(2, 2);
    Flocal_g = zeros(2, 1);

    n1 = elem_dofs(1);
    n2 = elem_dofs(2);

    x1 = node_coords(n1);
    x2 = node_coords(n2);

    u1 = soln_full(n1);
    u2 = soln_full(n2);
    u = [u1 u2];

    for gp = 1:nGP
        xi = gpts(gp);
        wt = gwts(gp);
        N = [0.5 * (1 - xi), 0.5 * (1 + xi)];
        dNdxi = [-0.5, 0.5];
        Jac = h / 2;
        dNdx = dNdxi / Jac;
        du = dNdx * u';
        x = N * [x1 x2]';
        % f = 10.0 * exp(-5 * x) - 4.0 * exp(-x); % subject to change
        f = 0;

        % advection
        Klocal = Klocal + (a * N' * dNdx) * Jac * wt;
        % reaction
        Klocal = Klocal + (s * N' * N) * Jac * wt;
        %diffusion
        Klocal = Klocal + (mu * dNdx' * dNdx) * Jac * wt;
        %force vector
        Flocal = Flocal + N' * f * Jac * wt;
        % Flocal = Flocal - N' * a * du * Jac * wt;
        % Flocal = Flocal - dNdx' * mu * du * Jac * wt;
    end

    Klocal_g = Klocal_g + Klocal;
    Flocal_g = Flocal_g + Flocal;
end

function Klocal_pg = petrovGalerkin(a, mu, h, alpha, nGP, gpts, gwts, elem_dofs, soln_full_pg, node_coords, Klocal, Flocal)
    Klocal_pg = zeros(2, 2);
    %     Flocal_pg = zeros(2, 1);

    n1 = elem_dofs(1);
    n2 = elem_dofs(2);

    x1 = node_coords(n1);
    x2 = node_coords(n2);

    u1 = soln_full_pg(n1);
    u2 = soln_full_pg(n2);
    u = [u1 u2];

    for gp = 1:nGP
        xi = gpts(gp);
        wt = gwts(gp);

        N = [0.5 * (1 - xi), 0.5 * (1 + xi)];
        dNdxi = [-0.5, 0.5];

        Jac = h / 2;

        dNdx = dNdxi / Jac;
        du = dNdx * u';

        x = N * [x1 x2]';
        f = 10.0 * exp(-5 * x) - 4.0 * exp(-x);
        % petrov galerkin advection
        Klocal = Klocal + (alpha * a * dNdx' * dNdx) * Jac * wt;
        % forcing term
        % Flocal = Flocal + (alpha * a * dNdx' * dNdx) * Jac * wt;
    end

    Klocal_pg = Klocal_pg + Klocal;
    %     Flocal_pg = Flocal_pg + Flocal;
end

function [Klocal_supg, Flocal_supg] = supg(a, mu, h, alpha, tau, s, nGP, gpts, gwts, elem_dofs, node_coords, soln_full_supg, Klocal, Flocal)
    Klocal_supg = zeros(2, 2);
    Flocal_supg = zeros(2, 1);

    n1 = elem_dofs(1);
    n2 = elem_dofs(2);

    x1 = node_coords(n1);
    x2 = node_coords(n2);

    u1 = soln_full_supg(n1);
    u2 = soln_full_supg(n2);
    u = [u1 u2];

    for gp = 1:nGP
        xi = gpts(gp);
        wt = gwts(gp);

        N = [0.5 * (1 - xi), 0.5 * (1 + xi)];
        dNdxi = [-0.5, 0.5];

        Jac = h / 2;

        dNdx = dNdxi / Jac;
        du = dNdx * u';
        uh = N * u';
        x = N * [x1 x2]';
        % f = 10.0 * exp(-5 * x) - 4.0 * exp(-x);
        f = 0;

        mod_test = a * dNdx' + abs(s) * N';

        % stabilization
        % Klocal = Klocal + tau * a ^ 2 * (dNdx' * dNdx) * Jac * wt;

        % reaction terms
        % Klocal = Klocal + tau * a^2 * (dNdx' * dNdx) * Jac * wt ...
        % + tau * a * s * (dNdx' * N) * Jac * wt;
        Klocal = Klocal + tau * (mod_test * (a * dNdx + s * N) * Jac * wt);

        % force vector
        res = a * du + s * uh - f;
        % Flocal = Flocal + tau * a * dNdx' * f * Jac * wt;
        Flocal = Flocal + tau * (mod_test * f) * Jac * wt;
        Flocal = Flocal - tau * dNdx' * a ^ 2 * du * Jac * wt;
    end

    Klocal_supg = Klocal_supg + Klocal;
    Flocal_supg = Flocal_supg + Flocal;

end

function [Klocal_ppv, Flocal_ppv] = ppv(a, mu, h, alpha, tau, s, nGP, gpts, gwts, elem_dofs, node_coords, soln_full_ppv, Klocal, Flocal)
    Klocal_ppv = zeros(2, 2);
    Flocal_ppv = zeros(2, 1);

    n1 = elem_dofs(1);
    n2 = elem_dofs(2);

    x1 = node_coords(n1);
    x2 = node_coords(n2);

    u1 = soln_full_ppv(n1);
    u2 = soln_full_ppv(n2);
    u = [u1 u2];

    eps = 1e-12;

    for gp = 1:nGP
        xi = gpts(gp);
        wt = gwts(gp);

        N = [0.5 * (1 - xi), 0.5 * (1 + xi)];
        dNdxi = [-0.5, 0.5];

        Jac = h / 2;

        dNdx = dNdxi / Jac;

        du = dNdx * u';
        uh = N * u';

        x = N * [x1 x2]';
        % f = 10.0 * exp(-5 * x) - 4.0 * exp(-x);
        f = 0;

        if (norm(du) < 1e-8)
            % res_ratio = (abs( s * uh - f) / (abs(du) + eps));
            res_ratio = abs((a * dNdx + s * N)) / abs((dNdx) + eps);

            chi = 2 / ((abs(s) * h) + (2 * abs(a)));

            kadd = max((abs(a - tau * a * s + tau * a * abs(s)) * h / 2) ...
                - (mu + tau * a ^ 2) + (s + tau * s * abs(s)) * h ^ 2/6, 0);

            % kadd = max((abs(a)* h/2) - (mu) + (s/6 * h^2), 0);
        else
            chi = 0;
            res_ratio = 0;
            kadd = 0;
        end

        % stabilization
        Klocal = Klocal + chi * res_ratio * kadd * (dNdx' * dNdx) * Jac * wt;

        Flocal = Flocal - res_ratio * chi * kadd * dNdx' * du * Jac * wt;
    end

    Klocal_ppv = Klocal_ppv + Klocal;
    Flocal_ppv = Flocal_ppv + Flocal;

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

function u_analytical = analyticalSolution(x, a, mu, s, L, uL, uR)

    f = 0;
    a = uL;
    b = uR;

    % jaiman paper
    m1 = (a + sqrt(a ^ 2 + 4 * mu * s)) / (2 * mu);

    m2 = (a - sqrt(a ^ 2 + 4 * mu * s)) / (2 * mu);

    num = a * exp(m2 * L) - ((f / s) * (exp(m2 * L) - 1)) - b;
    denum = exp(m2 * L) - exp(m1 * L);

    % c1 = (a * exp(m2 * l) - ((f / s) * (exp(m2 * l) - 1)) - b) / (exp(m2 * l) - exp(m1 * l));
    c1 = num / denum;
    c2 = a - c1 - (f / s);

    % u_analytical = (exp(a * x * mu) - 1) / (exp(pe) - 1);
    u_analytical = c1 * exp(m1 * x) + c2 * exp(m2 * x) + (f / s);

end
