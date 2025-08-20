function plot_results(U, gamma, elem_node_conn, node_coords)

q1 = U(1, :);
q2 = U(2, :) ./ q1;
q3 = U(3, :) ./ q1;
q4 = U(4, :);
u = [q2 q3];

p = zeros(max(size(U)), 1);

for i = 1:max(size(U))
    p(i) = (gamma - 1) * (q4(i) - 0.5 * q1(i) * norm([q2(i) q3(i)])^2);
end
% p = (gamma - 1) * (q4 - 0.5 * q1 .* norm(u).^2);

figure;
subplot(2, 2, 1);
% plot(x, rho, 'k-', 'LineWidth', 2);
% trimesh(elem_node_conn, node_coords(:, 1), node_coords(:, 2), rho);
trisurf(elem_node_conn, node_coords(:, 1), node_coords(:, 2), q1, "EdgeColor", "k", "FaceColor", "interp");
colormap jet; colorbar; view([0 0]);
title('Density');
xlabel('x'); zlabel('\rho');
grid on;
% ylim([0, 1.2]);

subplot(2, 2, 2);
% figure;
trisurf(elem_node_conn, node_coords(:, 1), node_coords(:, 2), q2, "EdgeColor", "k", "FaceColor", "interp");
colormap jet; colorbar; view([0 0]);
title('Velocity');
xlabel('x'); zlabel('x-velocity');
grid on;
%

subplot(2, 2, 3);
trisurf(elem_node_conn, node_coords(:, 1), node_coords(:, 2), p, "EdgeColor", "k", "FaceColor", "interp");
colormap jet; colorbar; view([0 0]);
title('Pressure');
xlabel('x'); ylabel('p');
grid on; 
% ylim([0, 1.2]);
%
subplot(2, 2, 4);
% figure;
trisurf(elem_node_conn, node_coords(:, 1), node_coords(:, 2), q4, "EdgeColor", "k", "FaceColor", "interp");
colormap jet; colorbar; view([0 0]);
title('Total Energy');
xlabel('x'); ylabel('E');
grid on;
%
% sgtitle('Riemann Shock Tube Solution (Taylor-Galerkin FEM)');
end


