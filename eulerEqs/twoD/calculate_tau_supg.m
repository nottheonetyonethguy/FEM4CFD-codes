function tau = calculate_tau_supg(U, gamma, dt, dNdx, dNdy, Jac)
q1 = U(1, :); % rho
q2 = U(2, :); % x momentum
q3 = U(3, :); % y momentum
q4 = U(4, :); % E

v1 = q2/q1; v2 = q3/q1;
v = [v1 v2];

p = (gamma - 1) * (q4 - 0.5 * q1 * norm(v)^2);
H = (q4 + p)/q1;
c = sqrt((gamma * p) / q1);

delN = [dNdx dNdy];

tau_sugn1 = (abs(v' * delN)) ^ (-1);

% tau_sugn1 = (c * norm(dNdx) + norm(U' * delN)) ^ (-1);
tau_sugn2 = dt / 2;

r = 2;
tau_ij = (1/tau_sugn1^r + 1/tau_sugn2^r) ^ (-1);

% g = 4; % gamma from spectral paper
% k = (1/Jac(1, 1))^g * (1/Jac(end,end))^g;
% k = 1;
tau = tau_ij;
% tau = k * tau_ij * eye(4);

end
