function [Klocal, Mlocal, K_supg, M_supg, K_shock] = stiffnessAndForce_quad4_euler2D(nodeNums, node_coords, elemData, matData, U, dt, soln_full, soln_full_prev)

% soln_full = u

npElem = size(nodeNums, 2); % number of points in an element

xNode = node_coords(nodeNums, 1); % x and y coordinates in the global frame
yNode = node_coords(nodeNums, 2); % x and y coordinates in the global frame

he = xNode(2) - xNode(1); % element length

gamma = matData(1, 1);

if (npElem == 3)
	nGP = 1;
	[gpts1, gpts2, gwts] = get_Gausspoints_tria(nGP);
	ELEMTYPE = 1;
elseif (npElem == 6)
	nGP = 3;
	[gpts1, gpts2, gwts] = get_Gausspoints_tria(nGP);
	ELEMTYPE = 1;
elseif (npElem == 4)
	nGP = 4;
	[gpts1, gpts2, gwts] = get_Gausspoints_quad(nGP);
	ELEMTYPE = 2;
elseif (npElem == 9)
	nGP = 9;
	[gpts1, gpts2, gwts] = get_Gausspoints_quad(nGP);
	ELEMTYPE = 2;
end

% discrete gradient-field matrix
G = zeros(2, npElem);

n_dof = 4; % no of dofs per node
n_nodes = 4; % no of nodes per element
matrix_size = n_nodes * n_dof;

Klocal = zeros(matrix_size, matrix_size);
Mlocal = zeros(matrix_size, matrix_size);

K_supg = zeros(matrix_size, matrix_size);
M_supg = zeros(matrix_size, matrix_size);

K_shock = zeros(matrix_size, matrix_size);

soln_elem = [U(:, 1); U(:, 2); U(:, 3); U(:, 4)];

param = [0.0; 0.0];

for gp = 1:nGP
	% get the coordinates and weights of quadrature point
	param(1) = gpts1(gp);
	param(2) = gpts2(gp);
	
	[N, dNdx, dNdy, Jac, detJ] = computeBasisFunctions2D(flag, ELEMTYPE, npElem, param, xNode, yNode);
	
	JacInv = inv(Jac);
	
	dvol = detJ * gwts(gp);
	
	% U_gp = U * N; % interpolated U values
	U_gp = zeros(n_dof, 1);
	for i = 1:n_nodes
		node_i_dofs = soln_elem((i-1)*n_dof + 1 : i*n_dof);
		U_gp = U_gp + N(i) * node_i_dofs;
	end
	
	tau = calculate_tau_supg(U_gp, gamma, dt, dNdx, dNdy, Jac);
	del = calculateDelShock(U_gp, gamma,dNdx, dNdy); 
    % del = 0;
	% tau = 1e-5;
	% tau = 0;
	[F1, F2, A1, A2] = calculateFluxJacobians(U_gp, gamma);
	
	for ii = 1:npElem
		G(1, ii) = dNdx(ii);
		G(2, ii) = dNdy(ii);
	end
	
	% mass matrix
	% Mlocal = (N * N') * eye(4) * dvol;
	
	% galerkin
	% Klocal_lhs = Klocal_lhs + N * (A1 * dNdx + A2 * dNdy)' * dvol;
	
	% assembly
	for i = 1:n_nodes
		for j = 1:n_nodes
			rowIdx = (i-1)*n_dof + 1 : i * n_dof;
			colIdx = (j-1)*n_dof + 1 : j * n_dof;
			
			K_ij = (N(i) * (A1 * dNdx(j) + A2 * dNdy(j))) * dvol; % K_local
			K_ij_supg = tau * (dNdx(i) * A1 + dNdy(i) * A2) * (A1 * dNdx(j) + A2 * dNdy(j)) * dvol;
			K_ij_shock = del * (dNdx(i) * dNdx(j) + dNdy(i) * dNdy(j)) * dvol * eye(4);
			
			M_ij = (N(i) * N(j) * eye(n_dof)) * dvol; % Mass matrix
			M_ij_supg = (tau * (A1 * dNdx(i) + A2 * dNdy(i))) * N(j) * dvol;
			
			Klocal(rowIdx, colIdx) = Klocal(rowIdx, colIdx) + K_ij;
			K_supg(rowIdx, colIdx) = K_supg(rowIdx, colIdx) + K_ij_supg;
			
			Mlocal(rowIdx, colIdx) = Mlocal(rowIdx, colIdx) + M_ij;
			M_supg(rowIdx, colIdx) = M_supg(rowIdx, colIdx) + M_ij_supg;
			
			K_shock(rowIdx, colIdx) = K_shock(rowIdx, colIdx) + K_ij_shock;
		end
	end
	
end
