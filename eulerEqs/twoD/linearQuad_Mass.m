function [Mlocal]=linearQuad_Mass(elmDat, IEN, e, XX, gp, gw)

degree = 1;
nlocal = 4;
nlbf   = 4;
ndof   = 2;
nsize  = nlocal*ndof;

sss        = elmDat(1);
finite     = (elmDat(2)>1);
bforce(1)  = elmDat(3);
bforce(2)  = elmDat(4);
rho        = elmDat(6);
E          = elmDat(7);
nu         = elmDat(8);

axsy=(sss==3);

thick = 1.0;

Mlocal=zeros(nsize,nsize); % Local mass matrix

xOrig=zeros(nlbf,1);
yOrig=zeros(nlbf,1);

xOrig(1) = XX(IEN(e,1),1);
xOrig(2) = XX(IEN(e,2),1);
xOrig(3) = XX(IEN(e,3),1);
xOrig(4) = XX(IEN(e,4),1);

yOrig(1) = XX(IEN(e,1),2);
yOrig(2) = XX(IEN(e,2),2);
yOrig(3) = XX(IEN(e,3),2);
yOrig(4) = XX(IEN(e,4),2);

nGP = max(size(gp));


for gp2=1:nGP

    JacMult = gw(gp2) * thick;
    param(2) = gp(gp2);
    
    for gp1=1:nGP

        param(1) = gp(gp1);

        [N, dN_dx, dN_dy, Jac] = computeBasisFunctions2D(0, 2, degree, param, xOrig, yOrig);
        xOrig
        yOrig

        fact = gw(gp1) * JacMult;

        dvol = Jac * fact;

        xx = 0.0;
        yy = 0.0;
        
        for ii=1:nlbf
            xx =  xx + N(ii)*xOrig(ii);
            yy =  yy + N(ii)*yOrig(ii);
        end

        if(axsy)
          dvol  = dvol  * 2.0*PI*yy;
        end

        for ii=1:nlbf
          bb3 = N(ii)*dvol*rho;

          TI   = 2*(ii-1)+1;
          TIp1 = TI+1;

          for jj=1:nlbf
              TJ   = 2*(jj-1)+1;
              TJp1 = TJ+1;

              Mlocal(TI,   TJ)     =  Mlocal(TI,   TJ)   + bb3*N(jj);
              Mlocal(TIp1, TJp1)   =  Mlocal(TIp1, TJp1) + bb3*N(jj);
          end
        end

    end%gp1
end%gp2





