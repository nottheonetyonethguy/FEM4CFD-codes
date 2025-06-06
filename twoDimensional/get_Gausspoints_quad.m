function [gpts1, gpts2, gwts] = get_Gausspoints_quad(nGP)

nGP1 = sqrt(nGP); % Number of points in each direction

% get quadrature points in each direction
[gpts_xi,  gwts_xi]   = get_Gausspoints_1D(nGP1);
[gpts_eta, gwts_eta] = get_Gausspoints_1D(nGP1);


gpts1 = zeros(nGP,1);
gpts2 = zeros(nGP,1);
gwts  = zeros(nGP,1);

gp = 1;
for gp2=1:nGP1
  for gp1=1:nGP1
    gpts2(gp) = gpts_eta(gp2);
    gpts1(gp) = gpts_xi(gp1);

    gwts(gp)  = gwts_eta(gp2)*gwts_xi(gp1);

    gp = gp+1;
  end
end

