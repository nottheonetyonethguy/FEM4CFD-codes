function [Klocal, Flocal] = calcStiffnessAndForce_1D2noded_AdvectionDiffusionReaction(nodeNums, node_coords, a, mu, fflag, f, soln_full, iteration, taufactor)

    n1 = nodeNums(1);
    n2 = nodeNums(2);

    x1 = node_coords(n1);
    x2 = node_coords(n2);

    u1 = soln_full(n1);
    u2 = soln_full(n2);

    nGP = 2; % what does nGP mean??

    [gpts, gwts] = get_Gausspoints_1D(nGP);

    Klocal = zeros(2, 2);
    Flocal = zeros(2, 1);

    N = zeros(2, 1);
    dNdxi = zeros(2, 1);
    dNdx = zeros(2, 1);

    he = x2 - x1;

    for gp = 1:nGP

	% get the xi coordinate and weight
	xi = gpts(gp);
	wt = gwts(gp);

	% compute the basis functions and their derivatives wrt xi
	N(1) = 0.5 * (1.0 - xi);
	N(2) = 0.5 * (1.0 + xi);

	dNdxi(1) = -0.5;
	dNdxi(2) = 0.5;

	% compute the Jacobian
	Jac = dNdxi(1) * x1 + dNdxi(2) * x2;

	% compute the derivatives of basis functions wrt to x
	dNdx(1) = dNdxi(1) / Jac;
	dNdx(2) = dNdxi(2) / Jac;

	% prepare G matrix
	G = [dNdx(1) dNdx(2)];

	xc = N(1) * x1 + N(2) * x2;

	u = N(1) * u1 + N(2) * u2;
	du = dNdx(1) * u1 + dNdx(2) * u2;

	if (fflag == 2)
	    a = xc;
	end

	% advection
	Klocal = Klocal + N * a * G * Jac * wt;

	% diffusion
	Klocal = Klocal + mu * G' * G * Jac * wt;

	% Example 1
	if (fflag == 1)
	    % f = 10.0 * exp(-5 * xc) - 4.0 * exp(-xc);
	    f = f;
	end

	% Hemker problem
	if (fflag == 2)
	    f = -mu * pi ^ 2 * cos(pi * xc) - pi * xc * sin(pi * xc);
	end

	res = a * du - f;

	% if(iteration > 1)
	    %   du
	    %   %res
	    % end

	Flocal = Flocal + N * f * Jac * wt;
	Flocal = Flocal - N * a * du * Jac * wt;
	Flocal = Flocal - G' * mu * du * Jac * wt;

	% Stabilisation
	%tau = 0.0*he/2/abs(a);
	Pe = a * he / (2.0 * mu);
	xi = 1.0 / tanh(Pe) - 1.0 / Pe;
	tau = xi * he / (2.0 * a) * taufactor;

	Klocal = Klocal + tau * a * G' * a * G * Jac * wt;

	Flocal = Flocal + tau * a * G' * f * Jac * wt;
	Flocal = Flocal - tau * a * G' * a * du * Jac * wt;

	% discontinuity capturing
	kadd = max(0.0, 0.5 * (a - tau * a * f + tau * a * abs(f)) * he - (mu + tau * a * a) + (f + tau * f * abs(f)) * he * he / 6.0);

	XXi = 2.0 / (abs(f) * he + 2 * abs(a));

	Rratio = abs(a * du - f) / abs(du);

	% ##  if(iteration > 1)
	    % ##    Klocal = Klocal + Rratio*XXi*kadd*G'*G*Jac*wt;
	    % ##    Flocal = Flocal - Rratio*XXi*kadd*G'*du*Jac*wt;
	    % ##  end

	aa = a;

	if (a > 0.0)
	    aa = -a;
	end
	% 
	Klocal = Klocal + aa * XXi * kadd * G' * G * Jac * wt;
	Flocal = Flocal - aa * XXi * kadd * G' * du * Jac * wt;

    end
