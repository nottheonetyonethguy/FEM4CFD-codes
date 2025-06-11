function  tau = calculate_stabilisation_term(a, b, mu)
% a and b are vectors
% mu is a scalar

anorm = norm(a);

if anorm == 0.0
	anorm = 1e-12;
end

bnorm = norm(b);

if bnorm == 0.0
	bnorm = 1e-12;
end


he = 2.0*anorm/bnorm;

Pe  =  anorm*he/(2.0*mu);
c1 = 1.0;
c2 = 1.0;
xi  =  c2/tanh(c1*Pe/c2) - c2/(c1*Pe);
tau =  xi/bnorm;


