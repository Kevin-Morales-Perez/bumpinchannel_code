function [uf] = rc_facev(pP,pA,dp,da,wlsq_op_P,wlsq_op_A,d_eps,eps_vec,normf_v,uP,uA)
%Rie - Chow interpolation face velocities
% This function calculates rie-chow face velocities (both components)
%using the nodal velocities and the gradient pressures
% inputs
%pP nodal values of pressures surrounding node P [pPW;pPN;pPE;pPS,pP];
%pA nodal values of pressures surrounding node A [pAW;pAN;pAE;pAS,pA];
%dp,da convection coeffitients at nodes P and A (Scalars)
%wlsq_op_P,wlsq_op_A weighted least squares operators at nodes P and A
%[2,4]
%normf_v norm unitary vector to face [1,2]
%eps unitary vector in the P_A direction [2,1]
%uP,uA nodal velocities at points P and A [u,v]
%d_eps distance between nodes (scalar)

%Output
%uf(Scalar)
%________________________computing pressure gradients______________________
pnodP=pP(5);
pnodA=pA(5);
delt_pP=pP(1:4)-pnodP;
delt_pA=pA(1:4)-pnodA;
gradpP=wlsq_op_P*delt_pP;
gradpA=wlsq_op_A*delt_pA;
eps_vec=eps_vec';
%________________________main formula _____________________________________
uf= 0.5*dot((uP+uA),normf_v) +0.5*(dp+da)*(pnodP - pnodA)/d_eps -0.5*dot(((dp*gradpP) + (da*gradpA)),eps_vec);
end
