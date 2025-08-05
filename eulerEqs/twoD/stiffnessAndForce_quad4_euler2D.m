 function [Klocal, Flocal] = stiffnessAndForce_quad4_AdvDiff2D(nodeNums, node_coords, force_applied, elemData, matData, soln_full, soln_full_prev)

 % soln_full = u

npElem = size(nodeNums, 2); % number of points in an element

xNode = node_coords(nodeNums, 1); % x and y coordinates in the global frame
yNode = node_coords(nodeNums, 2); % x and y coordinates in the global frame

he = xNode(2) - xNode(1); % element length

soln_elem = soln_full(nodeNums); % elemental solution, 4x1 vector
soln_elem_prev = soln_full_prev(nodeNums);

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

Klocal = zeros(8, 8);

param = [0.0; 0.0];

for gp = 1:nGP
	% get the coordinates and weights of quadrature point
	param(1) = gpts1(gp);
	param(2) = gpts2(gp);
	
	[N, dNdx, dNdy, Jac, detJ] = computeBasisFunctions2D(flag, ELEMTYPE, npElem, param, xNode, yNode);
	
	JacInv = inv(Jac);
	
	dvol = detJ * gwts(gp);
	
	for ii = 1:npElem
		G(1, ii) = dNdx(ii);
		G(2, ii) = dNdy(ii);
	end
	
	%du = zeros(2,1);
	du = G * soln_elem;
	du_prev = G * soln_elem_prev;
	
	% advection
	Klocal = Klocal + N * avec' * G * dvol;
	Flocal = Flocal - N * avec' * du * dvol;
	
	% reaction term
	Klocal = Klocal + s * G' * G * dvol;
	Flocal = Flocal - s * G' * du * dvol;
	
	% diffusion
	Klocal = Klocal + G' * mu * G * dvol;
	Flocal = Flocal - G' * mu * du * dvol;
	
	Flocal = Flocal + N * f * dvol;
	
	% Stabilisation
	% SUPG
	bvec = JacInv * avec;
	tau_e = calculate_stabilisation_term(avec, bvec, mu);
	% tau_e = calculate_tau_supg(avec, mu, s, xNode, yNode);
	
	% SUPG
	tau_SUPG = tau_e;
	Klocal = Klocal + tau_SUPG * G' * avec * avec' * G * dvol;
	
	% reaction term
	Klocal = Klocal + tau_SUPG * (avec' * G)' * (s * N)' * dvol;
	
	% combined GLS/SGS
	% mod_test = (avec' * G)' + abs(s) * N;
	% Klocal = Klocal + tau_SUPG * mod_test * (avec' * G + s * N') * dvol;
	%
	% reaction term
	% Klocal = Klocal + tau_SUPG * (avec' * G)' * (s * N)' * dvol;
	
	Flocal = Flocal - tau_SUPG * G' * avec * avec' * du * dvol;
	
	% Discontinuity capturing term
	% if(norm(du_prev) > 1.0e-6)
	%
	% 	% normalise du
	% 	du_prev = du_prev/norm(du_prev);
	%
	% 	aparllvec  = du_prev*dot(avec,du_prev);
	%
	% 	bparllvec = JacInv*aparllvec;
	%
	% 	tau_parll = calculate_stabilisation_term(aparllvec, bparllvec, mu);
	%
	% 	tau_DCO = max(0, tau_parll-tau_e);
	% 	%tau_DCO = tau_parll;
	% 	Klocal = Klocal + tau_DCO*G'*aparllvec*aparllvec'*G*dvol;
	% 	Flocal = Flocal - tau_DCO*G'*aparllvec*aparllvec'*du*dvol;
	% end
	
	% jaiman shock capturing
	if (norm(du_prev) > 1.0e-4)
		residual_vec = avec' * du_prev - mu * (dNdx' * dNdx + dNdy' * dNdy) * soln_elem_prev + s * N' * soln_elem_prev - f;
		
		v_mag = norm(avec); % advection velocity gradient
		grad_mag = norm(du_prev); % previous solution gradient
		
		if grad_mag < 1e-12
			grad_mag = 1e-12;
		end
		
		chi = 2 / (abs(s) * he + 2 * v_mag);
		
		ppv_multiplier = chi * abs(residual_vec') / grad_mag;
		
		% for streamline and crosswind directions
		if v_mag > 1e-12
			v_unit = avec / v_mag; % unit velocity
			P_streamline = (v_unit * v_unit'); %/ (v_unit' * v_unit); % v cross v / abs(v)^2
			P_crosswind = eye(2) - P_streamline;
			
			tau_ppv = calculate_tau_PPV(avec, mu, s, JacInv);
			
			k_add_s = max(abs(v_mag - tau_ppv * v_mag * s + tau_ppv * v_mag * abs(s)) * he / 2 - mu - tau_ppv * v_mag * v_mag + (s + tau_ppv * s * abs(s)) * (he ^ 2) / 6, 0);
			
			k_add_c = max(abs(v_mag + tau_ppv * v_mag * abs(s)) * he / 2 - mu + (s + tau_ppv * s * abs(s)) * (he ^ 2) / 6, 0);
			
			if k_add_s > 0.0
				K_ppv_s = ppv_multiplier * G' * k_add_s * P_streamline * G * dvol;
				F_ppv_s = ppv_multiplier * G' * k_add_s * P_streamline * du * dvol;
				
				Klocal = Klocal + K_ppv_s;
				Flocal = Flocal + F_ppv_s;
			end
			
			if k_add_c > 0.0
				K_ppv_c = ppv_multiplier * G' * k_add_c * P_crosswind * G * dvol;
				F_ppv_c = ppv_multiplier * G' * k_add_c * P_crosswind * du * dvol;
				
				Klocal = Klocal + K_ppv_c;
				Flocal = Flocal + F_ppv_c;
			end
			
		end
		
	end
	
	%Flocal = Flocal + tau1*avec*G'*f*dvol;
	
end
