function [del] = calculateDelShock(U, gamma, dNdx, dNdy)
q1 = U(1, :); % rho
q2 = U(2, :); % x momentum
q3 = U(3, :); % y momentum
q4 = U(4, :); % E

v1 = q2/q1; v2 = q3/q1;
p = (gamma - 1) * (q4 - 0.5 * q1 * (v1^2 + v2^2));
H = (q4 + p)/q1;
c = sqrt((gamma * p) / q1);

u = [v1 v2];
u_cha = c; u_int = u_cha;
rho_ref = 1;

h_ij = 0;
num = q1 * (dNdx + dNdy);
denum = norm(num);
j = num / denum;

for i = 1:4
		h_ij = h_ij + (c * abs(j(i) * (dNdx(i) + dNdy(i))));
end
h_shoc = 2 * 1 / h_ij;

tau_shoc = 0;
for i = 1:2
		beta = i;
		tau_shoc = tau_shoc + 0.5 * (h_shoc / (2 * u_cha)) * ((norm(q1 * (dNdx + dNdy)) * h_shoc) / rho_ref)^beta;
end
del = tau_shoc * u_int^2;
end
