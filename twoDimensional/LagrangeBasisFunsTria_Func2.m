function [N, dN_dxi1, dN_dxi2] = LagrangeBasisFunsTria_Func2(npElem, xi1, xi2)

    N = zeros(npElem,1);
    dN_dxi1 = zeros(npElem,1);
    dN_dxi2 = zeros(npElem,1);

    xi3 = 1.0 - xi1 - xi2;

    switch(npElem)
      case 1

          N(1) = 1.0;

          dN_dxi1(1) = 0.0;
          dN_dxi2(1) = 0.0;

      case 3

          N(1) = xi3;
          N(2) = xi1;
          N(3) = xi2;

          dN_dxi1(1) = -1.0;
          dN_dxi1(2) =  1.0;
          dN_dxi1(3) =  0.0;

          dN_dxi2(1) = -1.0;
          dN_dxi2(2) =  0.0;
          dN_dxi2(3) =  1.0;

      case 6

          N(1) = xi3*(2.0*xi3 - 1.0);
          N(2) = xi1*(2.0*xi1 - 1.0);
          N(3) = xi2*(2.0*xi2 - 1.0);
          N(4) = 4.0*xi1*xi3;
          N(5) = 4.0*xi1*xi2;
          N(6) = 4.0*xi2*xi3;

          dN_dxi1(1) = -4.0*xi3 + 1.0;
          dN_dxi1(2) =  4.0*xi1 - 1.0;
          dN_dxi1(3) =  0.0;
          dN_dxi1(4) =  4.0*(xi3 - xi1);
          dN_dxi1(5) =  4.0*xi2;
          dN_dxi1(6) = -4.0*xi2;

          dN_dxi2(1) = -4.0*xi3 + 1.0;
          dN_dxi2(2) =  0.0;
          dN_dxi2(3) =  4.0*xi2 - 1.0;
          dN_dxi2(4) = -4.0*xi1;
          dN_dxi2(5) =  4.0*xi1;
          dN_dxi2(6) =  4.0*(xi3 - xi2);

      otherwise

          printf("ERROR in LagrangeBasisFunsTria ... no basis functions defined for this npElem = %5d \n", npElem);

     end
