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
	e = (p / ((gamma - 1) * rho)) + 0.5 * u^2;
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
	K_supg_global = zeros(3*n_nodes, 3*n_nodes);
	K_sc_global = zeros(3*n_nodes, 3*n_nodes);
	M_supg_global = zeros(3*n_nodes, 3*n_nodes);
	% counter = 0;
	for elem = 1:n_elements
		node1 = elem; node2 = elem+1;
		UL = U(:, node1); UR = U(:, node2);
		[gpts, gwts] = get_Gausspoints_1D(nGP);
		Jac_elem = dx / 2;
		
		K_elem = zeros(6, 6);
		K_supg_elem = zeros(6, 6);
		K_sc_elem = zeros(6, 6);
		M_supg_elem = zeros(6, 6);
		% RHS_elem = zeros(6, 1);
		dof1 = (node1 - 1) * 3 + (1:3); dof2 = (node2 - 1) * 3 + (1:3);
		global_dofs = [dof1, dof2];
		
		for gp = 1:nGP
			xi = gpts(gp);
			wt = gwts(gp);
			
			N = [0.5 * (1 - xi), 0.5 * (1 + xi)];
			dN_dxi = [-0.5, 0.5];
			dN_dx = dN_dxi / Jac_elem;
			u_gp = N(1) * UL + N(2) * UR;
			% du_gp = [UL(2, :), UR(2, :)] * dN_dx';
			% stiffness matrix
			[F, A] = calculate_flux_jacobian(u_gp, gamma);
			tau = calculateTau_SUPG(UL, UR, dx, dt, Jac_elem, N, dN_dx, gamma);
			% tau = 0.002;
			for i = 1:2
				i_dofs = (i - 1) * 3 + (1:3);
				for j = 1:2
					j_dofs = (j - 1) * 3 + (1:3);
					K_ij = N(i) * A * dN_dx(j) * Jac_elem * wt;
					K_elem(i_dofs, j_dofs) = K_elem(i_dofs, j_dofs) + K_ij;
				end
			end
			% supg matrices
			for i = 1:2
				i_dofs = (i - 1) * 3 + (1:3);
				for j = 1:2
					j_dofs = (j - 1) * 3 + (1:3);
					M_ij_supg = tau * A * dN_dx(i) * N(j) * Jac_elem * wt;
					M_supg_elem(i_dofs, j_dofs) = M_supg_elem(i_dofs, j_dofs) + M_ij_supg;
					K_ij_supg = tau * (A * dN_dx(i))' * (A * dN_dx(j)) * Jac_elem * wt;
					K_supg_elem(i_dofs, j_dofs) = K_supg_elem(i_dofs, j_dofs) + K_ij_supg;
				end
			end
			% shock capturing
			for i = 1:2
				i_dofs = (i - 1) * 3 + (1:3);
				for j = 1:2
					j_dofs = (j - 1) * 3 + (1:3);
					K_ij_sc = tau * (dN_dx(i) * dN_dx(j)) * Jac_elem * wt * eye(3);
					K_sc_elem(i_dofs, j_dofs) = K_sc_elem(i_dofs, j_dofs) + K_ij_sc;
				end
			end
			
			% counter = counter + 1;
		end
		M_supg_global(global_dofs, global_dofs) = M_supg_global(global_dofs, global_dofs) + M_supg_elem;
		K_global(global_dofs, global_dofs) = K_global(global_dofs, global_dofs) + K_elem;
		K_supg_global(global_dofs, global_dofs) = K_supg_global(global_dofs, global_dofs) + K_supg_elem;
		K_sc_global(global_dofs, global_dofs) = K_sc_global(global_dofs, global_dofs) + K_sc_elem;
	end
	
	dofs_free = 4:(n_nodes * 3 - 3);
	
	LHS = (1/dt) * (M_global + M_supg_global) + (K_global + K_supg_global + K_sc_global);
	RHS = (1/dt) * (M_global + M_supg_global) * U(:);
	
	RHS = RHS - LHS(:, 1) * U(1, 1) - LHS(:, 2) * U(2, 1) - LHS(:, 3) * U(3, 1);
	RHS = RHS - LHS(:, end-2) * U(1, end) - LHS(:, end-1) * U(2, end) - LHS(:, end) * U(3, end);
	
	
	
	U_new = LHS(dofs_free, dofs_free) \ RHS(dofs_free);
	U(:, 2:end-1) = reshape(U_new, 3, n_nodes-2);
	
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


%% stabilization parameter tau
function tau = calculateTau_SUPG(UL, UR, dx, dt, Jac_elem, N, dN_dx, gamma)
% tau = zeros(3,3);

U = N(1) * UL + N(2) * UR; % u_gp
q1 = max(U(1, 1), 1e-10); % rho
q2 = U(2, 1); % momentum ,rho
q3 = U(3, 1); % energy, rho e


u = q2 / q1;
p = max((gamma - 1) * (q3 - 0.5 * q2 * u), 1e-10);
H = (q3 + p) / q1;

c = sqrt(gamma * p / q1);
h = dx;

du_dx = u * dN_dx(1) + u * dN_dx(2);

%% spectral
% tau_sugn1 = 1/(abs(u) * abs(dN_dx(1))) + (abs(u)* abs(dN_dx(2)));
tau_sugn2 = dt/2;
k = 1/Jac_elem;
if u == 0
	% tau_sugn1 = 0;
	tau_i = (1/tau_sugn2^2)^(-0.5);
else
	tau_sugn1 = 1/(norm(du_dx));
	tau_i = (1/tau_sugn1^2 + (1/tau_sugn2)^2)^(-0.5);
end

% tau_i = h / 2* (abs(u) + c);
tau = tau_i * eye(3);
% % tau = h / (abs(u) + c);

%% donea
% [F, A] = calculate_flux_jacobian(U, gamma);
% if u == 0
% 	lambda = diag([-c c]);
% 	R = [1, 1;...
% 		-c, c;...
% 		H, H];
% 	R_inv = inv(R' * R) * R';
% else
% 	lambda = diag([u - c u u + c]); % eigen vectors
% 	R = [1, 1, 1;...
% 		u-c, u, u+c;...
% 		H - u*c, 0.5 * u^2, H +  u*c];
% 	R_inv = inv(R);
% end
% a = R * abs(lambda) / R;
% [R, lambda] = eig(A);
% tau = (0.5 * h) * (R * abs(inv(lambda)) * inv(R));
% tau = (0.5 * h) * (R * inv(abs(lambda)) * R_inv);

%% catabriga

end

% function delta_hp = calculateSCterm(U, U_old, U_L, dx, dt, Jac_elem, dN_dx, gamma)
% end

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
