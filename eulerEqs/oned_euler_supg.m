% parameters
gamma = 1.4;
x_left = 0.0;
x_right = 1.0;
x_interface = 0.50;

n_elements = 100;
dx = (x_right - x_left) / n_elements;
n_nodes = n_elements + 1;

nGP = 2;

t_final = 0.2;
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
U_L = U(:, 1); U_R = U(:, end);

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
end
% M_lumped = diag(sum(M_global, 2));

% time stepping
t = 0.0;
dt = 0.001;
step = 0;
while t < t_final
	if t + dt > t_final
		dt = t_final - t;
	end
	
	K_global = zeros(3*n_nodes, 3*n_nodes);
	for elem = 1:n_elements
		node1 = elem; node2 = elem+1;
		UL = U(:, node1); UR = U(:, node2);
		[gpts, gwts] = get_Gausspoints_1D(nGP);
		Jac_elem = dx / 2;
		
		K_elem = zeros(6, 6); RHS_elem = zeros(6, 1);
		dof1 = (node1 - 1) * 3 + (1:3); dof2 = (node2 - 1) * 3 + (1:3);
		global_dofs = [dof1, dof2];
		
		for gp = 1:nGP
			xi = gpts(gp);
			wt = gwts(gp);
			
			N = [0.5 * (1 - xi), 0.5 * (1 + xi)];
			dN_dxi = [-0.5, 0.5];
			dN_dx = dN_dxi / Jac_elem;
			u_gp = N(1) * UL + N(2) * UR;
			% stiffness matrix
			[F, A] = calculate_flux_jacobian(u_gp, gamma);
			tau = calculateTau_SUPG(u_gp, dx, dt, gamma);
			for i = 1:2
				i_dofs = (i - 1) * 3 + (1:3);
				for j = 1:2
					j_dofs = (j - 1) * 3 + (1:3);
					K_ij = N(i) * A * dN_dx(j) * Jac_elem * wt;
					K_ij_supg = tau * (A * dN_dx(i))' * (A *dN_dx(j)) * Jac_elem * wt;
					K_elem(i_dofs, j_dofs) = K_elem(i_dofs, j_dofs) + K_ij + K_ij_supg;
				end
			end
		end
		residual = (u_gp - u_old_gp) / dt + A * (u_gp * dN_dx);
		K_global(global_dofs, global_dofs) = K_global(global_dofs, global_dofs) + K_elem;
	end
	
	LHS = (1/dt) * M_global + K_global;
	RHS = (1/dt) * M_global * U(:);
	
	U_new = LHS \ RHS;
	U = reshape(U_new, 3, n_nodes);
	
	U(: ,1) = U_L; U(:, end) = U_R;
	
	if any(~isfinite(U(:)))
		fprintf('NaN/Inf detected at step %d, time %.6f\n', step, t);
		break;
	end
	
	% time update
	t = t + dt;
	step = step + 1;
	fprintf('Step %d: Time = %.4f, dt = %.6f\n', step, t, dt);
end

plot_results(x, U, gamma);

%% Flux and Jacobian
function [F, A] = calculate_flux_jacobian(U, gamma)
F = zeros(3, 1);
A = zeros(3, 3);

q1 = max(U(1, 1), 1e-10); % rho
q2 = U(2, 1); % momentum ,rho u
q3 = U(3, 1); % energy, rho e

u = q2 / q1;
p = (gamma - 1) * (q3 - 0.5 * q2 * u);
H = (q3 + p) / q1;

F = [q2; q2 * u + p; (q3 + p) * u];
A = [0, 1, 0; ...
	0.5 * (gamma - 3) * u ^ 2, (3 - gamma) * u, gamma - 1; ...
	u * (0.5 * (gamma - 1) * u ^ 2 - H), H - (gamma - 1) * u ^ 2, gamma * u];

end

function R_elem = force_matrices(U_elem, dx, dt, gamma, nGP)
[gpts, gwts] = get_Gausspoints_1D(nGP);

R_elem = zeros(3, 3);
% K_global = zeros(3,3);
% K_supg = zeros(3,3);
Jac_elem = dx / 2;

for gp = 1:nGP
	xi = gpts(gp);
	wt = gwts(gp);
	
	N = [0.5 * (1 - xi), 0.5 * (1 + xi)];
	dN_dxi = [-0.5, 0.5];
	dN_dx = dN_dxi / Jac_elem;
	
	[F_elem, A_elem] = calculate_flux_jacobian(U_elem, gamma);
	% A_gp = 0.5 * (A_elem(:, :, 1) + A_elem(:, :, 2));
	A_gp = A_elem(:, :, gp);
	
	% taylor-galerkin terms:
	% stiffness matrix K
	K_global = A_gp * (N * dN_dx') * Jac_elem * wt;
	
	% supg terms
	tau = calculateTau_SUPG(U_elem, dx, dt, gamma);
	K_supg = tau * A_gp * ((N * dN_dx')/dt + A_gp * (dN_dx * dN_dx')) * Jac_elem * wt;
	% dF_dx = F_elem * dN_dx';
	% supg_term = tau * A_gp' * dN_dx' * dF_dx;
	% rhs = -dt * dF_dx + (dt^2 / 2) * A_gp * dF_dx + supg_term;
	R_elem = R_elem + K_global + K_supg;
	% assembly
	% for i = 1:2
	% 	i_dofs = (i-1)*3 + (1:3);
	%   R_elem(i_dofs) = R_elem(i_dofs) + rhs * N(i) * Jac_elem * wt;
	% end
end
end

%% stabilization parameter tau
function tau = calculateTau_SUPG(U, dx, dt, gamma)
tau = zeros(3,3);
rho = max(U(1, :), 1e-10);
mom = U(2, :);
E = max(U(3, :), 1e-10);

u = mom / rho;
p = (gamma - 1) * (E - 0.5 * rho * u^2);
p = max(p, 1e-10);
c = sqrt(gamma * p / rho);
h = dx;

tau_sugn1 = h ./(2*(c + (abs(u))));
tau_sugn2 = dt/2;
k = 2/h;
tau_i = k * (1/tau_sugn1^2 + 1/tau_sugn2^2)^(-0.5);
tau = tau_i * eye(3);
end

%% Plot results
function plot_results(x, U, gamma)
rho = U(1, :);
u = U(2, :) ./ rho;
E = U(3, :);
p = (gamma - 1) * (E - 0.5 * rho .* u.^2);

figure;
subplot(2, 2, 1);
plot(x, rho, 'k-', 'LineWidth', 2);
title('Density');
xlabel('x'); ylabel('\rho');
grid on; ylim([0, 1.2]);

subplot(2, 2, 2);
plot(x, u, 'k-', 'LineWidth', 2);
title('Velocity');
xlabel('x'); ylabel('u');
grid on;

subplot(2, 2, 3);
plot(x, p, 'k-', 'LineWidth', 2);
title('Pressure');
xlabel('x'); ylabel('p');
grid on; ylim([0, 1.2]);

subplot(2, 2, 4);
plot(x, E, 'm-', 'LineWidth', 2);
title('Total Energy');
xlabel('x'); ylabel('E');
grid on;

sgtitle('Riemann Shock Tube Solution (Taylor-Galerkin FEM)');
end
