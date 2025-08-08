function [F1, F2, A1, A2] = calculateFluxJacobian(U, gamma)
q1 = U(1, :); % rho
q2 = U(2, :); % x momentum
q3 = U(3, :); % y momentum
q4 = U(4, :); % E

u = q2 / q1; v = q3 / q1;
p = (gamma - 1) * (q4 - 0.5 * q1 * (u^2 + v^2));
H = q4 + p/q1;

F1 = [q2; q2 * u + p; q2 * v + p; H * q2];
F2 = [q3; q3 * u + p; q3 * v + p; H * q3];

A1 = [0, 1, 0, 0; ...
	0.5 * (gamma - 1) * (u^2 + v^2) - u^2, (3 - gamma) * u, (1 - gamma) * v, (gamma - 1);...
	- u * v, v, u, 0;...
	0.5 * (gamma - 1) * u * (u^2 + v^2) - u * H, - v^2 * (gamma - 1) + H, (1 - gamma) * u * v, gamma * u];

A2 = [0, 0, 1, 0; ...
	- u * v, v, u, 0;...
	0.5 * (gamma - 1) * (u^2 + v^2) - v^2, (1 - gamma) * u, (3 - gamma) * v, (gamma - 1);...
	0.5 * (gamma - 1) * v * (u^2 + v^2) - v * H, (1 - gamma) * u * v, - v^2 * (gamma - 1) + H, gamma * u];

end
