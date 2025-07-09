%% parameters
gamma = 1.4;
x_left = 0.0;
x_right = 1.0;
x_interface = 0.5;

n_elements = 100;
dx = (x_right - x_left) / n_elements;
n_nodes = n_elements + 1;

nGP = 2;
t_final = 0.2;

rho_L = 1.0; u_L = 0.0; p_L = 1.0;
rho_R = 0.125; u_R = 0.0; p_R = 0.1;

x = linspace(x_left, x_right, n_nodes);
U = zeros(3, n_nodes);
U_old = U;
% initial conditions
for i = 1:n_nodes
	if x(i) <= x_interface
		rho = rho_L; u_gp = u_L; p = p_L;
	else
		rho = rho_R; u_gp = u_R; p = p_R;
	end
	e = p / ((gamma - 1) * rho) + 0.5 * u_gp^2;
	U(:, i) = [rho; rho * u_gp; rho * e];
end
U_L = U(:, 1); U_R = U(:, end);

% mass matrix
M_global = zeros(3 * n_nodes, 3 * n_nodes);
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
M_lumped = diag(sum(M_global, 2));

% time stepping
t = 0.0;
dt = 0.0015;
step = 0;
while t < t_final
	if t + dt > t_final
		dt = t_final - t;
	end
	
	RHS = zeros(3 * n_nodes, 1);
	for elem = 1:n_elements
		UL = U(:, elem);
		UR = U(:, elem+1);
		
		[gpts, gwts] = get_Gausspoints_1D(nGP);
		for gp = 1:nGP
			xi = gpts(gp);
			wt = gwts(gp);
			
			N = [0.5 * (1 - xi), 0.5 * (1 + xi)];
			dN_dxi = [-0.5, 0.5];
			Jac_elem = dx / 2 ;
			dN_dx = dN_dxi / Jac_elem;
			u_gp = N(1) * UL + N(2) * UR;
			
			[F, A] = calculate_flux_jacobian(u_gp, gamma);
			% flux
            % F = N(1) * F + N(2) * F;
            % A = N(1) * A + N(2) * A;
			flux_term = (F * dN_dx) * Jac_elem * wt;
			dofs_1 = 3*(elem-1)+1:3*elem;
			dofs_2 = 3*elem+1:3*(elem+1);
			RHS(dofs_1) = RHS(dofs_1) + flux_term(:, 1);
			RHS(dofs_2) = RHS(dofs_2) + flux_term(:, 2);
			% second term
			ux = [UL, UR] * dN_dx';
			visc = (A^2) * ux;
			visc_term = 0.5 * dt * (visc * dN_dx) * Jac_elem * wt;
			RHS(dofs_1) = RHS(dofs_1) - visc_term(:, 1);
			RHS(dofs_2) = RHS(dofs_2) - visc_term(:, 2);
		end
	end
	% Ul = N(1) * U(:,1) + N(2) * U(:, 2);
    % Ur = N(1) * U(:, end) + N(2) * U(:, end);
    % % Ul = U(:, 1);
    % % Ur = U(:, end);
	% [Fl, Al] = calculate_flux_jacobian(Ul, gamma); % left boundary
	% [Fr, Ar] = calculate_flux_jacobian(Ur, gamma); % right boundary
	% % % boundary_term = br_t - bl_t;
	% RHS(1:3) = RHS(1:3) - (Fl + dt * 0.5 * (A * Fl));
	% RHS(end-2:end) = RHS(end-2:end) + (Fr + dt * 0.5 * (A * Fr));
	
	delta_U = M_lumped \ RHS;
	U = U + dt * reshape(delta_U, 3, n_nodes);
    U(:, 1) = U_L; U(:, end) = U_R;
	if any(~isfinite(U(:)))
		fprintf('NaN/Inf detected at step %d, time %.6f\n', step, t);
        break;
		U = U_old;
		dt = dt * 0.5;
		if dt < 1e-10
			error('Time step too small, simulation failed');
            break;
		end
		continue;
	end
	
	% U_old = U;
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
H = (q3 + p)/q1;

F = [q2; q2 * u + p; (q3 + p) * u];
A = [0, 1, 0; ...
	0.5 * (gamma - 3)* u^2, (3 - gamma) * u, gamma - 1; ...
	u * (0.5 * (gamma - 1) * u^2 - H), H - (gamma - 1) * u^2, gamma * u];

end

%% Plot results
function plot_results(x, U, gamma)
rho = U(1, :);
u = U(2, :) ./ rho;
E = U(3, :);
p = (gamma - 1) * (E - 0.5 * rho .* u.^2);

figure;
plot(x, rho, 'k-', 'LineWidth', 2);
title('Density');
xlabel('x'); ylabel('\rho');
grid on; ylim([0, 1.2]);

figure;
plot(x, u, 'k-', 'LineWidth', 2);
title('Velocity');
xlabel('x'); ylabel('u');
grid on;

% figure;
% plot(x, p, 'k-', 'LineWidth', 2);
% title('Pressure');
% xlabel('x'); ylabel('p');
% grid on; ylim([0, 1.2]);
%
% subplot(2, 2, 4);
figure;
plot(x, E, 'm-', 'LineWidth', 2);
title('Total Energy');
xlabel('x'); ylabel('E');
grid on;

% title('Riemann Shock Tube Solution (Taylor-Galerkin FEM)');
end
