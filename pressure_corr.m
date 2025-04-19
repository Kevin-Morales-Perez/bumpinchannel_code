function [p_star] = pressure_corr(p_m,len_f,d_eps_vec,eps_vec,normf_vec,u_vec,d_k,wl_op_mat)
%Pressure correction
%Pressure correction for collocated grid for SIMPLE algorithm
%p_m Matrix which containst the pressure nodal values needed for gradients
%p_m = 5-[0   0   pNN  0     0]
%      4-[0   pNW pN   PNE   0]
%      3-[pWW pW  pP   pE  pEE]
%      2-[0   pSW pS   pSE   0]
%      1-[0   0   pSS  0     0]
%        _____________________
%        1   2    3   4     5    

%pA=[pWest;pNorth;pEast;pSouth;Pprop]
%len_f Vector with length of faces [1,4]
%len_f =[len_fw,len_fn,len_fe,len_fs] 
%Least square operators, wlsq_op_W, wlsq_op_N,wlsq_op_E,wlsq_op_S,wlsq_op_P,
%[2,4]
%dp , da  coeffitients
%d_eps_v=[dW_eps,dN_eps,dE_eps,dS_eps]  distances between nodes

%eps_vec=[epsW_vec;epsW_vec;epsW_vec;epsW_vec]
%epsA_vec unit vector wich joins nodes P to A [1,2]

%normf_vec= [normfW_vec;normfN_vec;normfE_vec;normfS_vec]
%normfA_vec vector normal to face [1,2]

%u_vec =[u_W;u_N;u_E;u_S;u_P] 
%=u_A velocity at node A with components u_A=[uAx,vAy]
%d_k=[d_W,d_N,d_E,d_S] coeffitients

%wl_op_mat=[wl_op_w;wl_op_n;wl_op_e;wl_op_s;wl_op_p]
% weighted least squares operator matrix [10,4];
%wl_op_w weighted least operator [2,4]


pP=[p_m(3,2);p_m(4,3);p_m(3,4);p_m(2,3);p_m(3,3)];
pW=[p_m(3,1);p_m(4,2);p_m(3,3);p_m(2,2);p_m(3,2)];
pN=[p_m(4,3);p_m(5,3);p_m(4,4);p_m(3,3);p_m(4,3)];
pE=[p_m(3,3);p_m(4,4);p_m(3,5);p_m(2,4);p_m(3,4)];
pS=[p_m(2,2);p_m(3,3);p_m(2,4);p_m(1,3);p_m(2,3)];

wlsq_op_W=wl_op_mat(1:2,:);
wlsq_op_N=wl_op_mat(3:4,:);
wlsq_op_S=wl_op_mat(5:6,:);
wlsq_op_E=wl_op_mat(7:8,:);
wlsq_op_P=wl_op_mat(9:10,:);

dp=1;%not yet defined
da=1;%not yet defined

len_fw =len_f(1);
len_fn =len_f(2);
len_fe =len_f(3);
len_fs =len_f(4);

dW_eps=d_eps_vec(1);
dN_eps=d_eps_vec(2);
dE_eps=d_eps_vec(3);
dS_eps=d_eps_vec(4);

epsW_vec=eps_vec(1,:);
epsN_vec=eps_vec(2,:);
epsE_vec=eps_vec(3,:);
epsS_vec=eps_vec(4,:);

normfW_vec = normf_vec(1,:);
normfN_vec = normf_vec(2,:);
normfE_vec = normf_vec(3,:);
normfS_vec = normf_vec(4,:);

u_W=u_vec(1,:);
u_N=u_vec(2,:);
u_E=u_vec(3,:);
u_S=u_vec(4,:);
u_P=u_vec(5,:);

d_W =d_k(1);
d_N =d_k(2);
d_E =d_k(3);
d_S =d_k(4);



%___________west face______________
%normal velocity to face W.
%Rie - Chow interpolations in the internal part of the domain
uf_w = rc_facev(pP,pW,dp,da,wlsq_op_P,wlsq_op_W,dW_eps,epsW_vec,normfW_vec,u_P,u_W);
%Central diferencing for Nodes near to the boundaries
%mass flux at face W
massf_w=uf_w*len_fw;

%___________North face______________
%normal velocity to face N.
%Rie - Chow interpolations in the internal part of the domain
uf_n = rc_facev(pP,pN,dp,da,wlsq_op_P,wlsq_op_N,dN_eps,epsN_vec,normfN_vec,u_P,u_N);
%Central diferencing for Nodes near to the boundaries
%mass flux at face N
massf_n=uf_n*len_fn;

%___________East face______________
%normal velocity to face E.
%Rie - Chow interpolations in the internal part of the domain
uf_e = rc_facev(pP,pE,dp,da,wlsq_op_P,wlsq_op_E,dE_eps,epsE_vec,normfE_vec,u_P,u_E);
%Central diferencing for Nodes near to the boundaries
%mass flux at face E
massf_e=uf_e*len_fe;

%___________South face______________
%normal velocity to face S.
%Rie - Chow interpolations in the internal part of the domain
uf_s = rc_facev(pP,pS,dp,da,wlsq_op_P,wlsq_op_S,dS_eps,epsS_vec,normfS_vec,u_P,u_S);
%Central diferencing for Nodes near to the boundaries
%mass flux at face W
massf_s=uf_s*len_fs;

mass_inb=massf_w+massf_n+massf_e+massf_s;

%_____________Main formula___________________________________
p_W=p_m(3,2);
p_N=p_m(4,3);
p_E=p_m(3,4);
p_S=p_m(2,3);

a_p=len_fw*d_W + len_fn*d_N - len_fe*d_E -len_fs*d_S;

p_star= 1/a_p*(len_fw*d_W*p_W + len_fn*d_N*p_N - len_fe*d_E*p_E -len_fs*d_S*p_S-mass_inb);

end