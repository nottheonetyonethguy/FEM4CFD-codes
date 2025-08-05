function [gpts1, gpts2, gwts] = get_Gausspoints_tria(nGP)

gpts1 = zeros(nGP,1); % xi coordinate
gpts2 = zeros(nGP,1); % eta coordinate
gwts  = zeros(nGP,1); % weight

if(nGP == 1)
    gpts1(1) = 1.0/3.0;   gpts2(1) = 1.0/3.0;    gwts(1) = 0.5;
elseif(nGP == 3)
    gpts1(1) = 1.0/6.0;   gpts2(1) = 4.0/6.0;    gwts(1) = 1.0/6.0;
    gpts1(2) = 4.0/6.0;   gpts2(2) = 4.0/6.0;    gwts(2) = 1.0/6.0;
    gpts1(3) = 4.0/6.0;   gpts2(3) = 1.0/6.0;    gwts(3) = 1.0/6.0;
elseif(nGP == 4)
    gpts1(1) = 1.0/3.0;   gpts2(1) = 1.0/3.0;    gwts(1) = -27.0/96.0;
    gpts1(2) = 0.6;       gpts2(2) = 0.2;        gwts(2) =  25.0/96.0;
    gpts1(3) = 0.2;       gpts2(3) = 0.6;        gwts(3) =  25.0/96.0;
    gpts1(4) = 0.2;       gpts2(4) = 0.2;        gwts(4) =  25.0/96.0;
elseif(nGP == 7)
    gpts1(1) = 1.0/3.0;              gpts2(1) = 1.0/3.0;               gwts(1) =  0.1125;

    gpts1(2) = 0.797426985353;       gpts2(2) = 0.101286507323;        gwts(2) =  0.062969590272;
    gpts1(3) = 0.101286507323;       gpts2(3) = 0.797426985353;        gwts(3) =  0.062969590272;
    gpts1(4) = 0.101286507323;       gpts2(4) = 0.101286507323;        gwts(4) =  0.062969590272;

    gpts1(5) = 0.059715871789;       gpts2(5) = 0.470142064105;        gwts(5) =  0.066197076394;
    gpts1(6) = 0.470142064105;       gpts2(6) = 0.059715871789;        gwts(6) =  0.066197076394;
    gpts1(7) = 0.470142064105;       gpts2(7) = 0.470142064105;        gwts(7) =  0.066197076394;
end
