function tau = calculate_tau_supg(U, gamma, dt, dNdx, dNdy, Jac)
q1 = U(1, :); % rho
q2 = U(2, :); % x momentum
q3 = U(3, :); % y momentum
q4 = U(4, :); % E

v1 = q2/q1; v2 = q3/q1;
p = (gamma - 1) * (q4 - 0.5 * q1 * (v1^2 + v2^2));
H = (q4 + p)/q1;
c = sqrt((gamma * p) / q1);

tau_ij = 0;

num = q1 * (dNdx + dNdy);
denum = norm(num);
j = num / denum;

for i = 1:4
		tau_ij = tau_ij + abs(v1 * dNdx(i) + v2 * dNdy(2))... % first option
				+ c * abs(j(i) * (dNdx(i) + dNdy(i)))...% second option
				;
end

tau_sugn1 = 1/ tau_ij;
tau_sugn2 = dt / 2;
% tau_sugn1 = (c * norm(dNdx) + norm(U' * delN)) ^ (-1);

r = 2;
tau = (1/tau_sugn1^r + 1/tau_sugn2^r) ^ (-1/r);

% g = 4; % gamma from spectral paper
% k = (1/Jac(1, 1))^g * (1/Jac(end,end))^g;
% k = 1;
% tau = k * tau_ij * eye(4);

end
