function [N] = LagrangeBasisFunsQuad_Func1(npElem, xi1, xi2)

%   4---------3     4----7----3
%   |         |     |    |    |
%   |         |     |    |    |
%   |         |     8----9----6
%   |         |     |    |    |
%   |         |     |    |    |
%   1---------2     1----5----2
%
    N = zeros(npElem,1);

    switch(npElem)
        case 1

            N(1) = 1.0;

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

        case 9

            [Nu, dNu] = Lagrange_BasisFuns1D(2, xi1);
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

        otherwise

          printf("ERROR in LagrangeBasisFunsQuad... no basis functions defined for this npElem = %5d \n", npElem);

    end


