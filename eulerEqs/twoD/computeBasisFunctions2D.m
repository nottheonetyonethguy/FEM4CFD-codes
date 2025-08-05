function [N, dN_dx, dN_dy, Jac, detJ] = computeBasisFunctions2D(flag, ELEMTYPE, npElem, param, xNode, yNode)

N=zeros(npElem,1);
dN_dx=zeros(npElem,1);
dN_dy=zeros(npElem,1);

if(ELEMTYPE == 1) % triangular elements
    [N, dN_du1, dN_du2] = LagrangeBasisFunsTria_Func2(npElem, param(1), param(2));
else  % quad elements
    [N, dN_du1, dN_du2] = LagrangeBasisFunsQuad_Func2(npElem, param(1), param(2));
end

%Gradient of mapping from parameter space to physical space
Jac=zeros(2,2);
for ii=1:npElem
    Jac(1,1) = Jac(1,1) + xNode(ii) * dN_du1(ii) ;
    Jac(2,1) = Jac(2,1) + xNode(ii) * dN_du2(ii) ;
    Jac(1,2) = Jac(1,2) + yNode(ii) * dN_du1(ii) ;
    Jac(2,2) = Jac(2,2) + yNode(ii) * dN_du2(ii) ;
end

detJ  = det(Jac);

JacInv = inv(Jac);

%Compute derivatives of basis functions w.r.t physical coordinates
for ii=1:npElem
    dN_dx(ii) = dN_du1(ii) * JacInv(1,1) + dN_du2(ii) * JacInv(1,2);
    dN_dy(ii) = dN_du1(ii) * JacInv(2,1) + dN_du2(ii) * JacInv(2,2);
end

