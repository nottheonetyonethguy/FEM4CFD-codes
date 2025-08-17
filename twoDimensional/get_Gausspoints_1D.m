function [gpts, gwts] = get_Gausspoints_1D(nGP)

gpts = zeros(nGP,1);
gwts = zeros(nGP,1);

if(nGP == 1)
    gpts(1) = 0.0;                   gwts(1)  = 2.0;
elseif(nGP == 2)
    gpts(1) = -1.0/sqrt(3.0);        gwts(1)  = 1.0;
    gpts(2) =  1.0/sqrt(3.0);        gwts(2)  = 1.0;
elseif(nGP == 3)
    gpts(1) = -sqrt(0.6);            gwts(1)  = 5.0/9.0;
    gpts(2) =  0.0;                  gwts(2)  = 8.0/9.0;
    gpts(3) =  sqrt(0.6);            gwts(3)  = 5.0/9.0;
elseif(nGP == 4)
    gpts(1) = -0.861136;             gwts(1)  = 0.347855;
    gpts(2) = -0.339981;             gwts(2)  = 0.652145;
    gpts(3) =  0.339981;             gwts(3)  = 0.652145;
    gpts(4) =  0.861136;             gwts(4)  = 0.347855;
end
