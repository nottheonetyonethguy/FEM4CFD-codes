% 1D Riemann Shock Tube Problem using Euler Taylor-Galerkin Method
% Based on Sod's problem from the paper

% Problem parameters
gamma = 1.4;                % Ratio of specific heats
x_left = 0.0;              % Left boundary
x_right = 1.0;            % Right boundary
x_interface = 0.5;         % Initial interface position

% Initial conditions (Sod's problem)
% Left state (high pressure)
rho_L = 1.0;               % Density
u_L = 0.0;                 % Velocity
p_L = 1.0;                 % Pressure

% Right state (low pressure)
rho_R = 0.125;             % Density
u_R = 0.0;                 % Velocity
p_R = 0.1;                 % Pressure

% Numerical parameters
n_elements = 100;          % Number of elements
dx = (x_right - x_left) / n_elements;  % Element size
n_nodes = n_elements + 1;  % Number of nodes

% Time stepping
% CFL = 0.5;                 % CFL number
CFL = 1;                 % CFL number
t_final = 14.35;             % Final time
niter = 3;                 % Iterations for consistent mass matrix
c_lapidus = 1.0;           % Lapidus artificial viscosity constant

% Initialize solution
x = linspace(x_left, x_right, n_nodes);
U = zeros(3, n_nodes);     % [density; momentum; total energy]

% Set initial conditions
for i = 1:n_nodes
	if x(i) <= x_interface
		rho = rho_L;
		u = u_L;
		p = p_L;
	else
		rho = rho_R;
		u = u_R;
		p = p_R;
	end
	
	e = p / ((gamma - 1) * rho) + 0.5 * u^2;  % Specific total energy
	
	U(1, i) = rho;              % Density
	U(2, i) = rho * u;          % Momentum
	U(3, i) = rho * e;          % Total energy
end

% Time integration
t = 0.0;
time_step = 0;
dt = 0.205;

% Storage for results
U_history = [];
t_history = [];

fprintf('Starting simulation...\n');

while t < t_final
	% Calculate time step
	% dt = calculate_timestep(U, dx, gamma, CFL);
	
	% Ensure we don't overshoot final time
	if t + dt > t_final
		dt = t_final - t;
	end
	
	% Take one time step
	U = euler_taylor_galerkin_step(U, dx, dt, gamma, niter, c_lapidus);
	
	% Update time
	t = t + dt;
	time_step = time_step + 1;
	
	% Store results at certain intervals
	% if mod(time_step, 10) == 0 || t >= t_final
	U_history(:, :, end+1) = U;
	t_history(end+1) = t;
	fprintf('Time step %d, t = %.4f, dt = %.6f\n', time_step, t, dt);
	% end
end

% Plot results
plot_results(x, U, U_history, t_history, gamma);

fprintf('Simulation completed!\n');

function dt = calculate_timestep(U, dx, gamma, CFL)
% Calculate maximum allowable time step based on CFL condition

n_nodes = size(U, 2);
max_speed = 0;

for i = 1:n_nodes
	rho = U(1, i);
	u = U(2, i) / rho;
	E = U(3, i);
	p = (gamma - 1) * (E - 0.5 * rho * u^2);
	
	% Speed of sound
	c = sqrt(gamma * p / rho);
	
	% Maximum characteristic speed
	speed = abs(u) + c;
	max_speed = max(max_speed, speed);
end

dt = CFL * dx / max_speed;
end

function U_new = euler_taylor_galerkin_step(U, dx, dt, gamma, niter, c_lapidus)
% Single time step using Euler Taylor-Galerkin method

n_nodes = size(U, 2);
n_elements = n_nodes - 1;

% Calculate fluxes and Jacobians at current time
[F, A, dAdU] = calculate_flux_jacobian(U, gamma);

% Calculate source term derivatives (zero for Euler equations)
dSdU = zeros(3, 3, n_nodes);
S = zeros(3, n_nodes);

% Assemble element contributions
U_new = U;

for elem = 1:n_elements
	% Element nodes
	node1 = elem;
	node2 = elem + 1;
	
	% Element solution and derivatives
	U_elem = [U(:, node1), U(:, node2)];
	F_elem = [F(:, node1), F(:, node2)];
	A_elem = cat(3, A(:, :, node1), A(:, :, node2));
	
	% Calculate element matrices and vectors
	[M_elem, P_elem] = element_matrices(U_elem, F_elem, A_elem, dAdU(:,:,node1:node2), ...
		S(:,node1:node2), dSdU(:,:,node1:node2), dx, dt);
	
	% Apply consistent mass matrix iteration
	delta_U = zeros(3, 2);
	M_lumped = diag(diag(M_elem));  % Lumped mass matrix
	
	for iter = 1:niter
		if iter == 1
			delta_U_old = zeros(3, 2);
		else
			delta_U_old = delta_U;
		end
		
		% Solve lumped system
		rhs = P_elem - M_elem * delta_U_old(:);
		delta_U(:) = M_lumped \ rhs;
	end
	
	% Update nodal values
	U_new(:, node1) = U_new(:, node1) + delta_U(:, 1);
	U_new(:, node2) = U_new(:, node2) + delta_U(:, 2);
end

% Apply artificial viscosity (Lapidus smoothing)
if c_lapidus > 0
	U_new = apply_artificial_viscosity(U_new, dx, c_lapidus);
end

% Ensure positivity of density and pressure
U_new = enforce_positivity(U_new, gamma);
end

function [M_elem, P_elem] = element_matrices(U_elem, F_elem, A_elem, dAdU_elem, S_elem, dSdU_elem, dx, dt)
% Calculate element mass matrix and load vector

% Linear shape functions and derivatives
N = [0.5, 0.5; 0.5, 0.5];  % Shape functions at integration points
dN_dx = [-1/dx; 1/dx];      % Shape function derivatives

% Element mass matrix (consistent)
M_elem = (dx/6) * [2*eye(3), eye(3); eye(3), 2*eye(3)];

% Calculate load vector components
P_elem = zeros(6, 1);

for gp = 1:2  % Gauss points
	if gp == 1
		xi = -1/sqrt(3);
		w = 1.0;
	else
		xi = 1/sqrt(3);
		w = 1.0;
	end
	
	% Shape functions at Gauss point
	N1 = 0.5 * (1 - xi);
	N2 = 0.5 * (1 + xi);
	N_gp = [N1; N2];
	
	% Interpolate quantities at Gauss point
	U_gp = U_elem * N_gp;
	F_gp = F_elem * N_gp;
	S_gp = S_elem * N_gp;
	
	% Average Jacobian
	A_gp = 0.5 * (A_elem(:,:,1) + A_elem(:,:,2));
	dSdU_gp = 0.5 * (dSdU_elem(:,:,1) + dSdU_elem(:,:,2));
	
	% Calculate flux derivative
	dF_dx = F_elem * dN_dx;
	
	% First order terms
	term1 = dt * (S_gp - A_gp * (U_elem * dN_dx));
	
	% Second order terms
	term2 = (dt^2/2) * (dSdU_gp * (S_gp - A_gp * (U_elem * dN_dx)) - ...
		A_gp * dF_dx);
	
	% Add to load vector
	P_elem(1:3) = P_elem(1:3) + w * dx/2 * N1 * (term1 + term2);
	P_elem(4:6) = P_elem(4:6) + w * dx/2 * N2 * (term1 + term2);
end
end

function [F, A, dAdU] = calculate_flux_jacobian(U, gamma)
% Calculate flux vector, Jacobian matrix, and Jacobian derivatives

n_nodes = size(U, 2);
F = zeros(3, n_nodes);
A = zeros(3, 3, n_nodes);
dAdU = zeros(3, 3, n_nodes);

for i = 1:n_nodes
	rho = U(1, i);
	mom = U(2, i);
	E = U(3, i);
	
	u = mom / rho;
	p = (gamma - 1) * (E - 0.5 * rho * u^2);
	H = (E + p) / rho;  % Total enthalpy
	
	% Flux vector
	F(:, i) = [mom; rho*u^2 + p; mom*H];
	
	% Jacobian matrix A = dF/dU
	A(:, :, i) = [0, 1, 0;
		(gamma-3)*u^2/2, (3-gamma)*u, gamma-1;
		u*((gamma-1)*u^2 - H), H - (gamma-1)*u^2, gamma*u];
	
	% For this implementation, we'll approximate dAdU as zero
	% In a full implementation, you would compute the derivatives of A
	dAdU(:, :, i) = zeros(3, 3);
end
end

function U_smooth = apply_artificial_viscosity(U, dx, c_lapidus)
% Apply Lapidus artificial viscosity

n_nodes = size(U, 2);
U_smooth = U;

for i = 2:n_nodes-1
	% Calculate velocity at node i
	rho = U(1, i);
	u = U(2, i) / rho;
	
	% Calculate velocity gradient
	rho_left = U(1, i-1);
	rho_right = U(1, i+1);
	u_left = U(2, i-1) / rho_left;
	u_right = U(2, i+1) / rho_right;
	
	du_dx = (u_right - u_left) / (2 * dx);
	
	% Apply smoothing proportional to velocity gradient
	epsilon = c_lapidus * dx^2 * abs(du_dx);
	
	for var = 1:3
		U_smooth(var, i) = U(var, i) + epsilon * ...
			(U(var, i+1) - 2*U(var, i) + U(var, i-1)) / dx^2;
	end
end
end

function U_pos = enforce_positivity(U, gamma)
% Ensure density and pressure remain positive

U_pos = U;
n_nodes = size(U, 2);

for i = 1:n_nodes
	% Ensure positive density
	if U_pos(1, i) <= 0
		U_pos(1, i) = 1e-10;
	end
	
	% Ensure positive pressure
	rho = U_pos(1, i);
	mom = U_pos(2, i);
	E = U_pos(3, i);
	
	p = (gamma - 1) * (E - 0.5 * mom^2 / rho);
	
	if p <= 0
		% Adjust total energy to maintain positive pressure
		u = mom / rho;
		p_min = 1e-10;
		U_pos(3, i) = p_min / (gamma - 1) + 0.5 * rho * u^2;
	end
end
end

function plot_results(x, U_final, U_history, t_history, gamma)
% Plot the final results and animation

% Convert conservative to primitive variables
n_nodes = length(x);
rho = U_final(1, :);
u = U_final(2, :) ./ rho;
E = U_final(3, :);
p = (gamma - 1) * (E - 0.5 * rho .* u.^2);

% Create plots
figure('Position', [100, 100, 1200, 800]);

subplot(2, 2, 1);
plot(x, rho, 'b-', 'LineWidth', 2);
xlabel('x');
ylabel('Density');
title('Density Distribution');
grid on;

subplot(2, 2, 2);
plot(x, u, 'r-', 'LineWidth', 2);
xlabel('x');
ylabel('Velocity');
title('Velocity Distribution');
grid on;

subplot(2, 2, 3);
plot(x, p, 'g-', 'LineWidth', 2);
xlabel('x');
ylabel('Pressure');
title('Pressure Distribution');
grid on;

subplot(2, 2, 4);
plot(x, E, 'm-', 'LineWidth', 2);
xlabel('x');
ylabel('Total Energy');
title('Total Energy Distribution');
grid on;

sgtitle(sprintf('Riemann Shock Tube Solution at t = %.2f', t_history(end)));

% Print some statistics
fprintf('\nFinal Results:\n');
fprintf('Density range: [%.6f, %.6f]\n', min(rho), max(rho));
fprintf('Velocity range: [%.6f, %.6f]\n', min(u), max(u));
fprintf('Pressure range: [%.6f, %.6f]\n', min(p), max(p));
fprintf('Energy range: [%.6f, %.6f]\n', min(E), max(E));
end





