function [N, dN_dxi1, dN_dxi2] = LagrangeBasisFunsQuad_Func2(npElem, xi1, xi2)

%   4---------3     4----7----3
%   |         |     |    |    |
%   |         |     |    |    |
%   |         |     8----9----6
%   |         |     |    |    |
%   |         |     |    |    |
%   1---------2     1----5----2
%
    N = zeros(npElem,1);
    dN_dxi1 = zeros(npElem,1);
    dN_dxi2 = zeros(npElem,1);

    switch(npElem)
        case 1

            N(1) = 1.0;

            dN_dxi1(1) = 0.0;
            dN_dxi2(1) = 0.0;

        case 4

            [Nu, dNu] = Lagrange_BasisFuns1D(2, xi1);
            [Nv, dNv] = Lagrange_BasisFuns1D(2, xi2);

            %
            % v2          u1,v2  u2,v2
            % v1          u1,v1  u2,v1
            %   u1 u2

            N(1) = Nv(1)*Nu(1);
            N(2) = Nv(1)*Nu(2);
            N(4) = Nv(2)*Nu(1);
            N(3) = Nv(2)*Nu(2);

            dN_dxi1(1) = Nv(1)*dNu(1);
            dN_dxi1(2) = Nv(1)*dNu(2);
            dN_dxi1(4) = Nv(2)*dNu(1);
            dN_dxi1(3) = Nv(2)*dNu(2);

            dN_dxi2(1) = dNv(1)*Nu(1);
            dN_dxi2(2) = dNv(1)*Nu(2);
            dN_dxi2(4) = dNv(2)*Nu(1);
            dN_dxi2(3) = dNv(2)*Nu(2);

        case 9

            [Nu, dNu] = Lagrange_BasisFuns1D(3, xi1);
            [Nv, dNv] = Lagrange_BasisFuns1D(3, xi2);

            %
            % v2          u1,v2  u3,v2  u2,v2
            % v3          u1,v3  u3,v3  u2,v3
            % v1          u1,v1  u3,v1  u2,v1
            %   u1 u3 u2

                  N(1) = Nv(1)*Nu(1);
                  N(5) = Nv(1)*Nu(3);
                  N(2) = Nv(1)*Nu(2);

                  N(8) = Nv(3)*Nu(1);
                  N(9) = Nv(3)*Nu(3);
                  N(6) = Nv(3)*Nu(2);

                  N(4) = Nv(2)*Nu(1);
                  N(7) = Nv(2)*Nu(3);
                  N(3) = Nv(2)*Nu(2);

            dN_dxi1(1) = Nv(1)*dNu(1);
            dN_dxi1(5) = Nv(1)*dNu(3);
            dN_dxi1(2) = Nv(1)*dNu(2);

            dN_dxi1(8) = Nv(3)*dNu(1);
            dN_dxi1(9) = Nv(3)*dNu(3);
            dN_dxi1(6) = Nv(3)*dNu(2);

            dN_dxi1(4) = Nv(2)*dNu(1);
            dN_dxi1(7) = Nv(2)*dNu(3);
            dN_dxi1(3) = Nv(2)*dNu(2);
            %
            dN_dxi2(1) = dNv(1)*Nu(1);
            dN_dxi2(5) = dNv(1)*Nu(3);
            dN_dxi2(2) = dNv(1)*Nu(2);

            dN_dxi2(8) = dNv(3)*Nu(1);
            dN_dxi2(9) = dNv(3)*Nu(3);
            dN_dxi2(6) = dNv(3)*Nu(2);

            dN_dxi2(4) = dNv(2)*Nu(1);
            dN_dxi2(7) = dNv(2)*Nu(3);
            dN_dxi2(3) = dNv(2)*Nu(2);

        otherwise

          printf("ERROR in LagrangeBasisFunsQuad... no basis functions defined for this npElem = %5d \n", npElem);

    end


