function tau = calculate_tau_supg(u, k, s, xNode, yNode)
u_mag = norm(u);

Ae = polyarea(xNode, yNode); % area
h = sqrt(Ae); % length

tau = 1 / sqrt( (2*u_mag/h)^2 + (4*k/h^2)^2 + s^2 );
end
