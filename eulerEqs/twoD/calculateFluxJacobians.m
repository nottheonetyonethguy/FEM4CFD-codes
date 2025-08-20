function [F1, F2, A1, A2] = calculateFluxJacobian(U, gamma)
q1 = U(1, :); % rho
q2 = U(2, :); % x momentum
q3 = U(3, :); % y momentum
q4 = U(4, :); % rho E

v1 = q2 / q1; v2 = q3 / q1;
v = [v1 v2];
% V = norm(v);
p = (gamma - 1) * (q4 - 0.5 * q1 * norm(v)^2);
H = (q4 + p)/q1;

F1 = [q2; q2 * v1 + p; q2 * v2 + p; H * q2];
F2 = [q3; q3 * v1 + p; q3 * v2 + p; H * q3];

A1 = [0, 1, 0, 0; ...
		0.5 * (gamma - 1) * norm(v)^2 - v1^2, (3 - gamma) * v1, (1 - gamma) * v2, (gamma - 1);...
		- v1 * v2, v2, v1, 0;...
		0.5 * (gamma - 1) * v1 * norm(v)^2 - v1 * H, - (v1^2 * (gamma - 1) - H), (1 - gamma) * v1 * v2, gamma * v1];

A2 = [0, 0, 1, 0; ...
		-v1 * v2, v2, v1, 0;...
		0.5 * (gamma - 1) * (norm(v)^2) - v2^2, (1 - gamma) * v1, (3 - gamma) * v2, (gamma - 1);...
		0.5 * (gamma - 1) * v2 * (norm(v)^2) - v2 * H, (1 - gamma) * v1 * v2, - (v2^2 * (gamma - 1) - H), gamma * v1];

end
