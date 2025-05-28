clear all;
close all;

L = 1.0;
u = 1.0;
k = 0.05;
s = 600;
N = 50;
h = L / N;

phi_left = 8.0;
phi_right = 3.0;

x = linspace(0, L, N + 1)';

tau = 1 / sqrt((2 * u / h) ^ 2 + (4 * k / h ^ 2) ^ 2 + s ^ 2);

K = zeros(N + 1, N + 1);
F = zeros(N + 1, 1);
ke = zeros(N + 1, N + 1);
f = zeros(N + 1, 1);

for e = 1:N
    nodes = [e, e + 1];
    x_e = x(nodes);

    A_conv = u * [-1/2, 1/2; -1/2, 1/2];
    A_diff = k / h * [1, -1; -1, 1];
    A_react = s * h / 6 * [2, 1; 1, 2];

    L_op = u * (1 / h) * [-1, 1; -1, 1] + abs(s) * h / 6 * [2, 1; 1, 2];
    A_stab = tau * L_op * (u * (1 / h) * [-1, 1; -1, 1] + s * h / 6 * [2, 1; 1, 2]);

    k_add = max((abs(u) - tau * abs(u) * s + tau * abs(u) * abs(s)) * h / 2 - (k + tau * u ^ 2) + (s + tau * s * abs(s)) * h ^ 2/6, 0);
    chi = 2 / (abs(s) * h + 2 * abs(u));
    A_pos = chi * k_add / h * [1, -1; -1, 1];
    A_e = A_conv + A_diff + A_react + A_stab + A_pos;

    K(nodes, nodes) = K(nodes, nodes) + A_e;
end

K(1, :) = 0; K(1, 1) = 1; F(1) = phi_left;
K(end, :) = 0; K(end, end) = 1; F(end) = phi_right;

phi = K \ F;

figure;
plot(x, phi, 'b-o', 'LineWidth', 2);
xlabel('x');
ylabel('\phi(x)');
title('PPV Solution of 1D ADR Equation (u=1, k=0.01, s=60)');
grid on;
