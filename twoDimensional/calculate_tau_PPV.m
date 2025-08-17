function tau_ppv = calculate_tau_PPV(avec, mu, s, JacInv)

G_tensor = JacInv' * JacInv;

vGv = avec' * G_tensor * avec;
G_double_dot = trace(G_tensor * G_tensor);
C_I = 12.0;

denom = vGv + C_I * mu^2 * G_double_dot + s^2;

if denom < 1e-14
	tau_ppv = 0.0;
else
	tau_ppv = 1/sqrt(denom);
end
end
