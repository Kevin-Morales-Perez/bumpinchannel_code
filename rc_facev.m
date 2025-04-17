function [uf] = rc_facev(pP,pA,dp,da,wlsq_op_P,wlsq_op_A,eps,normf_v,uP,uA)
%Rie - Chow interpolation face velocities
% This function calculates rie-chow face velocities (both components)
%using the nodal velocities and the gradient pressures
%pP nodal values of pressures surrounding node P [pPW;pPN;pPE;pPS,pP];
%pA nodal values of pressures surrounding node A [pAW;pAN;pAE;pAS,pA];
%dp,da convection coeffitients at nodes P and A (Scalars)
%wlsq_op_P,wlsq_op_A weighted least squares operators at nodes P and A
%[2,4]
%normf_v norm unitary vector to face [1,2]
%eps unitary vector in the P_A direction [2,1]
%uP,uA nodal velocities at points P and A [u,v]

%________________________computing pressure gradients
pnodP=pP(5);
pnodA=pA(5);
delt_pP=pP(1:4)-pnodP;
delt_pA=pA(1:4)-pnodA;
gradpP=wlsq_op_P*delt_pP;
gradpA=wlsq_op_A*delt_pA;
%________________________main formula 
uf= 0.5*dot((uP+uA),normf_v) +0.5*(dp+da)*(pnodP + pnodA)/eps -0.5*dot(((dp*gradpP) + (da*gradpA)),eps);
end