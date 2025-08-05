function [N, dN_dxi]=Lagrange_BasisFuns1D(npElem, xi)

switch(npElem)
    case 1
        N(1) = 1.0;
        dN_dxi(1) = 0.0;

    case 2

        N(1) = 0.5*(1.0 - xi);
        N(2) = 0.5*(1.0 + xi);

        dN_dxi(1) = -0.5;
        dN_dxi(2) =  0.5;

    case 3

        val1 = xi*xi;

        N(1) = 0.5*(val1 - xi);
        N(2) = 0.5*(val1 + xi);
        N(3) = 1.0 - val1;

        val1 = 2.0*xi;

        dN_dxi(1) = 0.5*(val1 - 1.0);
        dN_dxi(2) = 0.5*(val1 + 1.0);
        dN_dxi(3) = -val1;

    otherwise

          printf('no basis functions defined for this degree = %5d \n', npElem);
end


