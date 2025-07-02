% parameters
gamma = 1.4;
x_left = 0.0;
x_right = 10.0;
x_interface = 5.0;

n_elements = 100;
dx = (x_right - x_left) / n_elements;
n_nodes = n_elements + 1;

nGP = 2;

% CFL = 0.5773;
t_final = 0.2;
% niter = 3;
c_lapidus = 0.1;

rho_L = 1.0; u_L = 0.0; p_L = 1.0;
rho_R = 0.125; u_R = 0.0; p_R = 0.1;

x = linspace(x_left, x_right, n_nodes);
U = zeros(3, n_nodes);

% initial conditions
for i = 1:n_nodes
	if x(i) <= x_interface
		rho = rho_L; u = u_L; p = p_L;
	else
		rho = rho_R; u = u_R; p = p_R;
	end
	e = p / ((gamma - 1) * rho) + 0.5 * u^2;
	U(:, i) = [rho; rho * u; rho * e];
end

% mass matrix
M_global = zeros(3*n_nodes, 3*n_nodes);
for elem = 1:n_elements
	[gpts, gwts] = get_Gausspoints_1D(nGP);
	M_elem = zeros(6, 6);
	
	Jac_elem = dx / 2;
	
	node1 = elem;
	node2 = elem + 1;
	dof1 = (node1-1)*3 + (1:3);
	dof2 = (node2-1)*3 + (1:3);
	global_dofs = [dof1, dof2];
	
	for gp = 1:nGP
		xi = gpts(gp);
		w = gwts(gp);
		
		N = [0.5 * (1 - xi), 0.5 * (1 + xi)];
		dN_dxi = [-0.5, 0.5];
		dN_dx = dN_dxi / Jac_elem;
		% assembly
		for i = 1:2
			i_dofs = (i-1)*3 + (1:3);
			for j = 1:2
				j_dofs = (j-1)*3 + (1:3);
				M_ij = N(i) * N(j) * eye(3) * Jac_elem * w;
				M_elem(i_dofs, j_dofs) = M_elem(i_dofs, j_dofs) + M_ij;
			end
		end
	end
	M_global(global_dofs, global_dofs) = M_global(global_dofs, global_dofs) + M_elem;
	M_lumped = diag(sum(M_global, 2));
end

% time stepping
t = 0.0;
dt = 0.001;
step = 0;
while t < t_final
	if t + dt > t_final
		dt = t_final - t;
	end
	
	U_old = U;
	U = euler_taylor_galerkin(U, dx, dt, gamma, c_lapidus, nGP, M_global, M_lumped);
	
	if any(~isfinite(U(:)))
		fprintf('NaN/Inf detected at step %d, time %.6f\n', step, t);
		U = U_old;
		dt = dt * 0.5;
		if dt < 1e-10
			error('Time step too small, simulation failed');
		end
		continue;
	end
	
	% time update
	t = t + dt;
	step = step + 1;
	fprintf('Step %d: Time = %.4f, dt = %.6f\n', step, t, dt);
end

plot_results(x, U, gamma);

%% calculate time step based on CFL condition
function dt = calculate_timestep(U, dx, gamma, CFL)
n_nodes = size(U, 2);
max_speed = 0;
for i = 1:n_nodes
	rho = U(1, i);
	if rho <= 0
		rho = 1e-10;
	end
	u = U(2, i) / rho;
	E = U(3, i);
	p = (gamma - 1) * (E - 0.5 * rho * u^2);
	if p <= 0
		p = 1e-10;
	end
	c = sqrt(gamma * p / rho);
	max_speed = max(max_speed, abs(u) + c);
end
dt = CFL * dx / max_speed;
end

%% euler taylor-galerkin method
function U_new = euler_taylor_galerkin(U, dx, dt, gamma, c_lapidus, nGP, M_global, M_lumped)
n_nodes = size(U, 2);
n_elements = n_nodes - 1;

% flux, jacobians
% [F, A] = calculate_flux_jacobian(U, gamma);

F_global = zeros(3*n_nodes, 1);

for elem = 1:n_elements
	node1 = elem;
	node2 = elem + 1;
	
	U_elem = U(:, node1:node2);
	% F_elem = F(:, node1:node2);
	% A_elem = cat(3, A(:, :, node1), A(:, :, node2));
	
	F_elem = force_matrices(U_elem, dx, dt, gamma, nGP);
	
	dof1 = (node1-1)*3 + (1:3);
	dof2 = (node2-1)*3 + (1:3);
	global_dofs = [dof1, dof2];
	
	F_global(global_dofs) = F_global(global_dofs) + F_elem;
end

delta_U_global = zeros(3*n_nodes, 1);

residual = F_global - M_global * delta_U_global;
delta_U_global = M_lumped \ residual;

U_new = U;
for i = 1:n_nodes
	dofs = (i-1)*3 + (1:3);
	U_new(:, i) = U(:, i) + delta_U_global(dofs);
end

U_new = apply_artificial_viscosity(U_new, dx, c_lapidus);

U_new = enforce_positivity(U_new, gamma);
end

function R_elem = force_matrices(U_elem, dx, dt, gamma, nGP)
[gpts, gwts] = get_Gausspoints_1D(nGP);

R_elem = zeros(6, 1);

Jac_elem = dx / 2;

for gp = 1:nGP
	xi = gpts(gp);
	w = gwts(gp);
	
	N = [0.5 * (1 - xi), 0.5 * (1 + xi)];
	dN_dxi = [-0.5, 0.5];
	dN_dx = dN_dxi / Jac_elem;
	
	% U_gp = U_elem * N';
	% F_gp = F_elem * N';
	% calc flux and jac inside the gauss loop, mass matrix outside the time
	% loop
	
	[F_elem, A_elem] = calculate_flux_jacobian(U_elem, gamma);
	% A_gp = 0.5 * (A_elem(:, :, 1) + A_elem(:, :, 2));
	A_gp = A_elem(:, :, gp);
	
	% taylor-galerkin terms:
	
	dF_dx = F_elem * dN_dx;
	
	rhs = -dt * dF_dx + (dt^2 / 2) * A_gp * dF_dx;
	
	% assembly
	for i = 1:2
		i_dofs = (i-1)*3 + (1:3);
		R_elem(i_dofs) = R_elem(i_dofs) + rhs * N(i) * Jac_elem * w;
	end
end
end

%% Flux and Jacobian
function [F, A] = calculate_flux_jacobian(U, gamma)
% n_nodes = size(U, 2);
n_nodes = 1;
F = zeros(3, n_nodes);
A = zeros(3, 3, n_nodes);

% for i = 1:n_nodes
i = 1;
rho = max(U(1, i), 1e-10);
mom = U(2, i);
E = max(U(3, i), 1e-10);

u = mom / rho;
p = (gamma - 1) * (E - 0.5 * rho * u^2);
p = max(p, 1e-10);
H = (E + p) / rho;

F(:, 1) = [mom; rho * u^2 + p; mom * H];

A(:, :, 1) = [0, 1, 0;
	0.5*(gamma-3)*u^2, (3-gamma)*u, gamma-1;
	u*(0.5*(gamma-1)*u^2 - H), H - (gamma-1)*u^2, gamma*u];
% end
end

%% artificial viscosity lapidus
function U_smooth = apply_artificial_viscosity(U, dx, c_lapidus)
n_nodes = size(U, 2);
U_smooth = U;

for i = 2:n_nodes-1
	rho_L = max(U(1, i-1), 1e-10);
	rho_R = max(U(1, i+1), 1e-10);
	rho_C = max(U(1, i), 1e-10);
	
	gamma = 1.4;
	p_L = (gamma - 1) * (U(3, i-1) - 0.5 * U(2, i-1)^2 / rho_L);
	p_R = (gamma - 1) * (U(3, i+1) - 0.5 * U(2, i+1)^2 / rho_R);
	p_C = (gamma - 1) * (U(3, i) - 0.5 * U(2, i)^2 / rho_C);
	
	dp_dx = abs(p_R - p_L) / (2 * dx);
	p_avg = (p_L + p_R + p_C) / 3;
	
	if dp_dx > 0.1 * p_avg / dx
		epsilon = c_lapidus * dx^2 * dp_dx;
		for var = 1:3
			U_smooth(var, i) = U(var, i) + epsilon * ...
				(U(var, i+1) - 2*U(var, i) + U(var, i-1)) / dx^2;
		end
	end
end
end

%% positivity (density and pressure)
function U_pos = enforce_positivity(U, gamma)
U_pos = U;
n_nodes = size(U, 2);

for i = 1:n_nodes
	rho_min = 1e-10;
	p_min = 1e-10;
	
	U_pos(1, i) = max(U_pos(1, i), rho_min);
	
	rho = U_pos(1, i);
	mom = U_pos(2, i);
	E = U_pos(3, i);
	p = (gamma - 1) * (E - 0.5 * mom^2 / rho);
	
	if p <= p_min
		u = mom / rho;
		U_pos(3, i) = p_min / (gamma - 1) + 0.5 * rho * u^2;
	end
end
end

%% Plot results
function plot_results(x, U, gamma)
rho = U(1, :);
u = U(2, :) ./ rho;
E = U(3, :);
p = (gamma - 1) * (E - 0.5 * rho .* u.^2);

figure;
% subplot(2, 2, 1);
plot(x, rho, 'k-', 'LineWidth', 2);
title('Density');
xlabel('x'); ylabel('\rho');
grid on; ylim([0, 1.2]);

% subplot(2, 2, 2);
figure;
plot(x, u, 'k-', 'LineWidth', 2);
title('Velocity');
xlabel('x'); ylabel('u');
grid on;

figure;
% subplot(2, 2, 3);
plot(x, p, 'k-', 'LineWidth', 2);
title('Pressure');
xlabel('x'); ylabel('p');
grid on; ylim([0, 1.2]);
%
% subplot(2, 2, 4);
% plot(x, E, 'm-', 'LineWidth', 2);
% title('Total Energy');
% xlabel('x'); ylabel('E');
% grid on;

% title('Riemann Shock Tube Solution (Taylor-Galerkin FEM)');
end
