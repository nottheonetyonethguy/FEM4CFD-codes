function [] = oneD_Stabilization(a, mu, h)
    %% petrov galerkin
    Pe = (a * h) / (2 * mu); % peclet number
    beta = coth(Pe) - (1 / Pe); % beta
    mudash = beta * (a * h) / 2;

    %% stream upwind
    %% supg
    %% stabilization
end
