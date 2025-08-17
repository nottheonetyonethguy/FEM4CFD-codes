function plot_results(x, U, gamma, elem_node_conn, node_coords)

rho = U(1, :);
u = U(2, :) ./ rho;
E = U(3, :);
p = (gamma - 1) * (E - 0.5 * rho .* u.^2);

figure;
% subplot(2, 2, 1);
% plot(x, rho, 'k-', 'LineWidth', 2);
% trimesh(elem_node_conn, node_coords(:, 1), node_coords(:, 2), rho);
trisurf(elem_node_conn, node_coords(:, 1), node_coords(:, 2), rho, "EdgeColor", "k", "FaceColor", "interp");
colormap jet; colorbar; view([0 0]);
title('Density');
xlabel('x'); zlabel('\rho');
grid on;
% ylim([0, 1.2]);

% subplot(2, 2, 2);
% plot(x, u, 'k-', 'LineWidth', 2);
% title('Velocity');
% xlabel('x'); ylabel('u');
% grid on;
%
% subplot(2, 2, 3);
% plot(x, p, 'k-', 'LineWidth', 2);
% title('Pressure');
% xlabel('x'); ylabel('p');
% grid on; ylim([0, 1.2]);
%
% subplot(2, 2, 4);
% plot(x, E, 'm-', 'LineWidth', 2);
% title('Total Energy');
% xlabel('x'); ylabel('E');
% grid on;
%
% sgtitle('Riemann Shock Tube Solution (Taylor-Galerkin FEM)');
end


